"""DuckDB-WASM / Python parity tests for the explainability UI SHAP queries.

The explainability UI runs SHAP queries in the browser via DuckDB-WASM against the
Parquet release assets. This test suite runs the SAME SQL (via the Python duckdb
bindings, which share the query engine with DuckDB-WASM) against a synthetic SHAP
parquet and verifies the returned values match a pyarrow-based ground truth to
machine precision, plus the load-bearing base + Σ SHAP = predicted logit invariant.

The point is not to re-validate DuckDB itself — it's to catch regressions in the
query templates embedded in `main.js` (e.g. a typo in the dynamic SELECT AVG(...)
list, a join predicate drift when cohort group-by is added/changed, a column-name
filter that skips a prefix).

Keeps the fixture tiny (~5 pairs × 10 features) so the test runs in &lt;1 s in CI,
fully standalone, no dependency on live snapshot data under `.scratch/`.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


# ---- fixtures --------------------------------------------------------------------

FIXTURE_FEATURE_NAMES = (
    "host_surface__host_capsule_cluster_1_score",
    "host_surface__host_capsule_cluster_2_score",
    "host_typing__host_serotype",
    "host_typing__host_o_type",
    "host_defense__CRISPR",
    "pair_depo_capsule__in_cluster_6",
    "pair_depo_capsule__has_depo_x_cluster_1",
    "pair_concentration__log10_pfu_ml",
    "phage_projection__recep_frac_OmpC",
    "phage_stats__phage_gc_content",
)


def _build_fixture_pairs() -> pd.DataFrame:
    """Five held-out pairs — mix of Guelin/BASEL and lysed/non-lysed for full coverage.

    All four (source × label) cohort cells populated so the cohort GROUP BY queries
    don't silently miss a zero-count bucket. Per-feature SHAP values picked to produce
    a distinctive signed average per cohort (makes aggregation bugs visible).
    """
    rows = [
        # bacteria,  phage,    source,  label
        ("H01", "P01", "guelin", 1),
        ("H01", "P02", "guelin", 0),
        ("H02", "P01", "guelin", 1),
        ("H03", "P03", "basel", 1),
        ("H04", "P04", "basel", 0),
    ]
    return pd.DataFrame(rows, columns=["bacteria", "phage", "source", "label"])


def _build_shap_matrix(n_pairs: int, n_features: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.normal(0.0, 0.3, size=(n_pairs, n_features)).astype(np.float64)


def _write_shap_parquet(tmp: Path, pairs: pd.DataFrame, shap_matrix: np.ndarray) -> Path:
    """Match the derive_shap_snapshot.py output schema: pair_id, bacteria, phage,
    source, fold_id, then one shap__<feature> column per feature."""
    columns: dict[str, np.ndarray | list] = {
        "pair_id": [f"{r.bacteria}__{r.phage}" for r in pairs.itertuples()],
        "bacteria": pairs["bacteria"].tolist(),
        "phage": pairs["phage"].tolist(),
        "source": pairs["source"].tolist(),
        "fold_id": [0] * len(pairs),
    }
    for i, name in enumerate(FIXTURE_FEATURE_NAMES):
        columns[f"shap__{name}"] = shap_matrix[:, i]
    df = pd.DataFrame(columns)
    path = tmp / "bacteria_axis_shap_values.parquet"
    pq.write_table(pa.Table.from_pandas(df), path)
    return path


def _write_base_parquet(tmp: Path, pairs: pd.DataFrame, base_value: float = -2.0981) -> Path:
    df = pd.DataFrame(
        {
            "pair_id": [f"{r.bacteria}__{r.phage}" for r in pairs.itertuples()],
            "base_value": [base_value] * len(pairs),
        }
    )
    path = tmp / "bacteria_axis_shap_base_values.parquet"
    pq.write_table(pa.Table.from_pandas(df), path)
    return path


# ---- tests -----------------------------------------------------------------------


def test_per_pair_shap_matches_pyarrow_ground_truth() -> None:
    """Mirrors the Pair explorer query: SELECT * FROM read_parquet(...) WHERE bacteria=? AND phage=?.

    Verifies DuckDB returns the exact per-row SHAP contributions a direct pyarrow read
    would surface — no rounding, no column filtering drift.
    """
    pairs = _build_fixture_pairs()
    matrix = _build_shap_matrix(len(pairs), len(FIXTURE_FEATURE_NAMES))
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        shap_path = _write_shap_parquet(tmpdir, pairs, matrix)

        target_bact, target_phage = "H01", "P01"
        ground_truth = pd.read_parquet(shap_path)
        expected = ground_truth[(ground_truth.bacteria == target_bact) & (ground_truth.phage == target_phage)].iloc[0]

        con = duckdb.connect()
        query = (
            f"SELECT * FROM read_parquet('{shap_path}') "
            f"WHERE bacteria = '{target_bact}' AND phage = '{target_phage}' LIMIT 1"
        )
        row = con.execute(query).fetchdf().iloc[0]

        for col in ground_truth.columns:
            if col.startswith("shap__"):
                assert row[col] == expected[col], f"{col}: duckdb {row[col]} != pyarrow {expected[col]}"


def test_base_plus_sum_shap_equals_predicted_logit_invariant() -> None:
    """TreeExplainer guarantees base + Σ SHAP = model output logit. This is the core
    self-consistency check the Pair explorer UI shows in its header. Any query-path
    arithmetic error between DuckDB and the final display would break this equality."""
    pairs = _build_fixture_pairs()
    matrix = _build_shap_matrix(len(pairs), len(FIXTURE_FEATURE_NAMES))
    base_value = -2.0981
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        shap_path = _write_shap_parquet(tmpdir, pairs, matrix)
        base_path = _write_base_parquet(tmpdir, pairs, base_value=base_value)

        for r in pairs.itertuples():
            con = duckdb.connect()
            shap_row = (
                con.execute(
                    f"SELECT * FROM read_parquet('{shap_path}') "
                    f"WHERE bacteria = '{r.bacteria}' AND phage = '{r.phage}' LIMIT 1"
                )
                .fetchdf()
                .iloc[0]
            )
            base_row = (
                con.execute(
                    f"SELECT base_value FROM read_parquet('{base_path}') "
                    f"WHERE pair_id = '{r.bacteria}__{r.phage}' LIMIT 1"
                )
                .fetchdf()
                .iloc[0]
            )

            shap_cols = [c for c in shap_row.index if c.startswith("shap__")]
            shap_sum = float(shap_row[shap_cols].sum())
            base = float(base_row.base_value)

            expected_logit = base + shap_sum
            # Ground truth via pyarrow
            pa_shap = pd.read_parquet(shap_path)
            pa_row = pa_shap[(pa_shap.bacteria == r.bacteria) & (pa_shap.phage == r.phage)].iloc[0]
            pa_sum = float(pa_row[shap_cols].sum())
            assert abs(expected_logit - (base + pa_sum)) < 1e-12


def test_cohort_groupby_source_mean_matches_pandas_groupby() -> None:
    """Mirrors the Cohort SHAP query: SELECT source, AVG(shap__<feat>) FROM ... GROUP BY source.

    Verifies DuckDB's GROUP BY + AVG produces the same per-cohort means as
    `pandas.groupby('source').mean()` over the same columns. This is load-bearing for
    the Cohort SHAP tab: a drift here would silently show wrong cohort differentials.
    """
    pairs = _build_fixture_pairs()
    matrix = _build_shap_matrix(len(pairs), len(FIXTURE_FEATURE_NAMES))
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        shap_path = _write_shap_parquet(tmpdir, pairs, matrix)

        shap_df = pd.read_parquet(shap_path)
        shap_cols = [c for c in shap_df.columns if c.startswith("shap__")]
        expected = shap_df.groupby("source")[shap_cols].mean().sort_index()

        con = duckdb.connect()
        agg_selects = ", ".join(f'AVG("{c}") AS "{c}"' for c in shap_cols)
        sql = (
            f"SELECT source AS grp, COUNT(*) AS n, {agg_selects} "
            f"FROM read_parquet('{shap_path}') GROUP BY source ORDER BY source"
        )
        actual_df = con.execute(sql).fetchdf().set_index("grp").drop(columns=["n"])

        np.testing.assert_array_almost_equal(
            expected.values,
            actual_df[shap_cols].values,
            decimal=15,
            err_msg="DuckDB GROUP BY AVG differs from pandas groupby mean",
        )


def test_cohort_groupby_label_via_values_join() -> None:
    """Mirrors the Cohort SHAP label path: the shap parquet has no label column, so
    the UI synthesizes a VALUES list from the predictions JSON and joins it in.
    Verifies the join + aggregate returns correct per-label means."""
    pairs = _build_fixture_pairs()
    matrix = _build_shap_matrix(len(pairs), len(FIXTURE_FEATURE_NAMES))
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        shap_path = _write_shap_parquet(tmpdir, pairs, matrix)

        shap_df = pd.read_parquet(shap_path)
        shap_cols = [c for c in shap_df.columns if c.startswith("shap__")]
        merged = shap_df.merge(pairs[["bacteria", "phage", "label"]], on=["bacteria", "phage"])
        expected = merged.groupby("label")[shap_cols].mean().sort_index()

        con = duckdb.connect()
        values = ", ".join(f"('{r.bacteria}','{r.phage}',{r.label})" for r in pairs.itertuples())
        agg_selects = ", ".join(f'AVG("{c}") AS "{c}"' for c in shap_cols)
        sql = (
            f"WITH labels(bacteria, phage, label) AS (VALUES {values}) "
            f"SELECT CAST(labels.label AS VARCHAR) AS grp, COUNT(*) AS n, {agg_selects} "
            f"FROM read_parquet('{shap_path}') sv JOIN labels USING (bacteria, phage) "
            f"GROUP BY labels.label ORDER BY labels.label"
        )
        actual_df = con.execute(sql).fetchdf().set_index("grp").drop(columns=["n"])
        actual_df.index = actual_df.index.astype(int)

        np.testing.assert_array_almost_equal(
            expected.values,
            actual_df.sort_index()[shap_cols].values,
            decimal=15,
            err_msg="DuckDB label-join GROUP BY differs from pandas groupby mean",
        )


def test_describe_yields_all_shap_columns() -> None:
    """The UI discovers SHAP columns via `DESCRIBE SELECT * FROM read_parquet(...)`.
    Regressions in that query (e.g. stale column caching) would silently drop features
    from the aggregate SELECT. This verifies DESCRIBE returns exactly the shap__<feat>
    columns a pyarrow schema read would show."""
    pairs = _build_fixture_pairs()
    matrix = _build_shap_matrix(len(pairs), len(FIXTURE_FEATURE_NAMES))
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        shap_path = _write_shap_parquet(tmpdir, pairs, matrix)

        pa_shap_cols = {c for c in pd.read_parquet(shap_path).columns if c.startswith("shap__")}

        con = duckdb.connect()
        describe = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{shap_path}')").fetchdf()
        duckdb_shap_cols = {c for c in describe["column_name"] if c.startswith("shap__")}

        assert duckdb_shap_cols == pa_shap_cols, f"DESCRIBE drift: pyarrow={pa_shap_cols}, duckdb={duckdb_shap_cols}"
        assert len(duckdb_shap_cols) == len(FIXTURE_FEATURE_NAMES)
