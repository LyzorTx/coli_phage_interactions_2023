from __future__ import annotations

import csv
import json
from pathlib import Path

import joblib
import pytest

from lyzortx.pipeline.autoresearch import build_contract, prepare_cache


def write_raw_interactions(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "bacteria",
        "bacteria_index",
        "phage",
        "image",
        "replicate",
        "plate",
        "log_dilution",
        "X",
        "Y",
        "score",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=";")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_fasta(path: Path, stem: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{stem}\nATGC\n", encoding="utf-8")


def write_json_file(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_locked_comparator_artifacts(root: Path) -> dict[str, object]:
    benchmark_summary_path = root / "track_g" / "tg02_gbm_calibration" / "tg02_benchmark_summary.json"
    feature_lock_path = root / "track_g" / "tg05_feature_subset_sweep" / "tg05_locked_v1_feature_config.json"
    model_summary_path = root / "track_g" / "tg01_v1_binary_classifier" / "tg01_model_summary.json"

    write_json_file(benchmark_summary_path, {"metric": "auprc", "value": 0.5})
    write_json_file(feature_lock_path, {"feature_blocks": ["adsorption"]})
    write_json_file(model_summary_path, {"model_type": "gbm"})

    return {
        "artifact_id": "track_g_clean_v1_locked_benchmark",
        "benchmark_summary_path": str(benchmark_summary_path),
        "feature_lock_path": str(feature_lock_path),
        "model_summary_path": str(model_summary_path),
        "evaluation_protocol_id": "steel_thread_v0_st03_split_v1",
        "locked_feature_blocks": ["adsorption"],
        "selection_source": "test comparator fixture",
    }


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_fixture_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for bacteria in ("B1", "B2", "B3", "B4", "B5"):
        for phage in ("P1", "P2"):
            for replicate, dilution, score in (
                ("1", "0", "1" if (bacteria, phage) in {("B1", "P1"), ("B2", "P2")} else "0"),
                ("2", "0", "0"),
                ("3", "0", "0"),
                ("1", "-1", "0"),
                ("2", "-1", "0"),
                ("1", "-2", "0"),
                ("2", "-2", "0"),
                ("3", "-2", "0"),
                ("1", "-4", "0"),
            ):
                rows.append(
                    {
                        "bacteria": bacteria,
                        "bacteria_index": bacteria,
                        "phage": phage,
                        "image": "",
                        "replicate": replicate,
                        "plate": "plate1",
                        "log_dilution": dilution,
                        "X": "0",
                        "Y": "0",
                        "score": score,
                    }
                )
    return rows


def test_build_top_level_schema_manifest_freezes_reserved_slots() -> None:
    schema = prepare_cache.build_top_level_schema_manifest()

    assert schema["supported_search_splits"] == ["train", "inner_val"]
    assert schema["disallowed_search_splits"] == ["holdout"]
    assert schema["slot_order"] == [
        "host_defense",
        "host_surface",
        "host_typing",
        "host_stats",
        "phage_projection",
        "phage_stats",
        "phage_kmer",
    ]
    assert schema["feature_slots"]["host_surface"] == {
        "entity_key": "bacteria",
        "join_keys": ["bacteria"],
        "column_family_prefix": "host_surface__",
        "block_role": "host",
        "reserved_feature_columns": [],
        "reserved_feature_column_count": 0,
    }
    assert schema["feature_slots"]["phage_projection"]["join_keys"] == ["phage"]


def test_main_writes_search_cache_without_holdout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    raw_path = tmp_path / "raw_interactions.csv"
    host_dir = tmp_path / "hosts"
    phage_dir = tmp_path / "phages"
    output_root = tmp_path / "outputs"
    cache_dir = output_root / "search_cache_v1"
    comparator_root = tmp_path / "comparator"
    tl17_output_dir = tmp_path / "tl17_runtime"

    write_raw_interactions(raw_path, build_fixture_rows())
    for bacteria in ("B1", "B2", "B3", "B4", "B5"):
        write_fasta(host_dir / f"{bacteria}.fna", bacteria)
    for phage in ("P1", "P2"):
        write_fasta(phage_dir / f"{phage}.fna", phage)
    tl17_output_dir.mkdir(parents=True)
    (tl17_output_dir / "tl17_rbp_reference_bank.faa").write_text(">ref_11\nMPEPTIDE\n", encoding="utf-8")
    (tl17_output_dir / "tl17_rbp_reference_metadata.csv").write_text(
        "reference_id,phage,family_id,gene_name,protein_index,annotation,phrog,protein_length_aa\n"
        "ref_11,P0,RBP_PHROG_11,P0_CDS_0001,1,tail fiber,11,8\n",
        encoding="utf-8",
    )
    (tl17_output_dir / "tl17_rbp_family_metadata.csv").write_text(
        "family_id,column_name,supporting_phage_count,supporting_reference_count\n"
        "RBP_PHROG_11,tl17_phage_rbp_family_11_percent_identity,2,1\n",
        encoding="utf-8",
    )
    write_json_file(
        tl17_output_dir / "schema_manifest.json",
        {
            "feature_block": "tl17_rbp_family_projection",
            "columns": [
                {"name": "phage", "dtype": "string"},
                {"name": "tl17_phage_rbp_family_11_percent_identity", "dtype": "float64"},
                {"name": "tl17_rbp_reference_hit_count", "dtype": "int64"},
            ],
            "family_score_columns": ["tl17_phage_rbp_family_11_percent_identity"],
            "reference_hit_count_column": "tl17_rbp_reference_hit_count",
        },
    )
    joblib.dump(
        {
            "family_rows": [
                {
                    "family_id": "RBP_PHROG_11",
                    "column_name": "tl17_phage_rbp_family_11_percent_identity",
                    "supporting_phage_count": 2,
                    "supporting_reference_count": 1,
                }
            ],
            "reference_rows": [
                {
                    "reference_id": "ref_11",
                    "phage": "P0",
                    "family_id": "RBP_PHROG_11",
                    "gene_name": "P0_CDS_0001",
                    "protein_index": 1,
                    "annotation": "tail fiber",
                    "phrog": "11",
                    "protein_sequence": "MPEPTIDE",
                }
            ],
            "matching_policy": {
                "min_query_coverage": 0.7,
                "mmseqs_command": ["mmseqs"],
            },
        },
        tl17_output_dir / "tl17_rbp_runtime.joblib",
    )
    write_json_file(
        tl17_output_dir / "tl17_phage_compatibility_manifest.json",
        {
            "task_id": "TL17",
            "matching_policy": {
                "min_query_coverage": 0.7,
                "mmseqs_command": ["mmseqs"],
            },
            "counts": {
                "retained_family_count": 1,
                "retained_reference_protein_count": 1,
            },
            "outputs": {
                "projected_feature_csv": str(tl17_output_dir / "tl17_panel_projected_phage_features.csv"),
            },
        },
    )

    monkeypatch.setattr(
        build_contract,
        "COMPARATOR_BENCHMARK",
        write_locked_comparator_artifacts(comparator_root),
    )
    fast_path_calls: list[dict[str, object]] = []

    def fake_build_host_surface_rows_fast_path(
        *,
        assemblies,
        output_dir,
        max_workers,
        include_lps_core_type,
    ) -> dict[str, object]:
        output_dir.mkdir(parents=True, exist_ok=True)
        fast_path_calls.append(
            {
                "assemblies": tuple(path.stem for path in assemblies),
                "output_dir": output_dir,
                "max_workers": max_workers,
                "include_lps_core_type": include_lps_core_type,
            }
        )
        return {
            "rows": [
                {
                    "bacteria": path.stem,
                    "host_o_antigen_type": f"{path.stem}_O",
                    "host_o_antigen_score": 10.0 + index,
                    "host_receptor_btub_score": float(index),
                    "host_capsule_profile_kpsc_score": 0.0,
                }
                for index, path in enumerate(sorted(assemblies), start=1)
            ],
            "schema": {
                "includes_lps_core_type": False,
            },
            "runtime_metadata": {
                "runtime_id": "test_fast_path",
                "legacy_nhmmer_path_forbidden": True,
                "total_elapsed_seconds": 1.0,
            },
        }

    monkeypatch.setattr(
        prepare_cache.run_all_host_surface,
        "build_host_surface_rows_fast_path",
        fake_build_host_surface_rows_fast_path,
    )

    def fake_build_host_defense_slot_artifact(
        *,
        args: object,
        cache_dir: Path,
        split_rows: dict[str, list[dict[str, str]]],
    ) -> dict[str, object]:
        slot_dir = cache_dir / "feature_slots" / "host_defense"
        slot_dir.mkdir(parents=True, exist_ok=True)
        retained_bacteria = sorted(
            {
                str(row["bacteria"])
                for rows in split_rows.values()
                for row in rows
                if str(row["retained_for_autoresearch"]) == "1"
            }
        )
        artifact_path = slot_dir / prepare_cache.SLOT_FEATURES_FILENAME
        with artifact_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["bacteria", "host_defense__AbiD"])
            writer.writeheader()
            for bacteria in retained_bacteria:
                writer.writerow({"bacteria": bacteria, "host_defense__AbiD": "1"})
        return {
            "artifact_path": str(artifact_path),
            "artifact_sha256": build_contract.sha256_file(artifact_path),
            "columns": ["host_defense__AbiD"],
            "column_count": 1,
            "entity_count": len(retained_bacteria),
            "build_manifest_path": str(slot_dir / prepare_cache.HOST_DEFENSE_BUILD_MANIFEST_FILENAME),
        }

    monkeypatch.setattr(
        prepare_cache,
        "_build_host_defense_slot_artifact",
        fake_build_host_defense_slot_artifact,
    )

    def fake_derive_host_typing_features(
        assembly_path: Path,
        *,
        bacteria_id: str | None = None,
        output_dir: Path,
        picard_metadata_path: Path | None,
    ) -> dict[str, object]:
        bacteria = bacteria_id or assembly_path.stem
        output_dir.mkdir(parents=True, exist_ok=True)
        feature_row = {
            "bacteria": bacteria,
            "host_clermont_phylo": f"{bacteria}_phylo",
            "host_st_warwick": f"{len(bacteria)}",
            "host_o_type": f"{bacteria}_O",
            "host_h_type": f"{bacteria}_H",
            "host_serotype": f"{bacteria}_O:{bacteria}_H",
        }
        return {
            "schema": {"feature_block": "host_typing"},
            "feature_row": feature_row,
            "comparison": None,
            "manifest": {
                "runtime_caveats": [
                    {
                        "bacteria": bacteria,
                        "caller": "sequence_type",
                        "field": "host_st_warwick",
                        "raw_value": "",
                        "normalized_value": "",
                        "message": "fixture caveat",
                    }
                ]
            },
        }

    monkeypatch.setattr(
        prepare_cache.derive_host_typing_features,
        "derive_host_typing_features",
        fake_derive_host_typing_features,
    )
    monkeypatch.setattr(
        prepare_cache.derive_host_typing_features,
        "build_host_typing_schema",
        lambda: {
            "caller_envs": {
                "phylogroup": "phylogroup_caller",
                "serotype": "serotype_caller",
                "sequence_type": "sequence_type_caller",
            }
        },
    )

    def fake_derive_host_stats_features(
        assembly_path: Path,
        *,
        bacteria_id: str | None = None,
        output_dir: Path,
    ) -> dict[str, object]:
        bacteria = bacteria_id or assembly_path.stem
        output_dir.mkdir(parents=True, exist_ok=True)
        return {
            "schema": {"feature_block": "host_stats"},
            "feature_row": {
                "bacteria": bacteria,
                "host_sequence_record_count": len(bacteria),
                "host_genome_length_nt": len(bacteria) * 100,
                "host_gc_content": round(len(bacteria) / 10, 3),
                "host_n50_contig_length_nt": len(bacteria) * 50,
            },
            "manifest": {},
        }

    monkeypatch.setattr(
        prepare_cache.derive_host_stats_features,
        "derive_host_stats_features",
        fake_derive_host_stats_features,
    )

    phage_projection_calls: list[dict[str, object]] = []

    def fake_project_phage_feature_rows_batched(
        phage_paths: list[Path],
        *,
        runtime_payload: dict[str, object],
        reference_fasta_path: Path,
        scratch_root: Path,
    ) -> list[dict[str, object]]:
        phage_projection_calls.append(
            {
                "phages": tuple(path.stem for path in phage_paths),
                "reference_fasta_path": reference_fasta_path,
                "scratch_root": scratch_root,
                "runtime_payload": runtime_payload,
            }
        )
        scratch_root.mkdir(parents=True, exist_ok=True)
        return [
            {
                "phage": path.stem,
                "tl17_phage_rbp_family_11_percent_identity": 70.0 + index,
                "tl17_rbp_reference_hit_count": index,
            }
            for index, path in enumerate(phage_paths, start=1)
        ]

    monkeypatch.setattr(
        prepare_cache.tl17_runtime,
        "project_phage_feature_rows_batched",
        fake_project_phage_feature_rows_batched,
    )
    monkeypatch.setattr(
        prepare_cache.tl17_runtime,
        "project_phage_feature_row",
        lambda *_args, **_kwargs: pytest.fail("Per-phage TL17 projection path should not be used"),
    )

    phage_stats_row_calls: list[dict[str, object]] = []

    def fake_build_phage_stats_feature_row(
        phage_fasta_path: Path,
        *,
        phage_id: str | None = None,
    ) -> dict[str, object]:
        phage = phage_id or phage_fasta_path.stem
        phage_stats_row_calls.append({"phage": phage, "path": phage_fasta_path})
        return {
            "phage": phage,
            "phage_sequence_record_count": len(phage),
            "phage_genome_length_nt": len(phage) * 100,
            "phage_gc_content": round(len(phage) / 10, 3),
            "phage_n50_contig_length_nt": len(phage) * 50,
        }

    monkeypatch.setattr(
        prepare_cache.derive_phage_stats_features,
        "build_phage_stats_feature_row",
        fake_build_phage_stats_feature_row,
    )
    monkeypatch.setattr(
        prepare_cache.derive_phage_stats_features,
        "derive_phage_stats_features",
        lambda *_args, **_kwargs: pytest.fail("Per-phage phage-stats artifact writes should not be used"),
    )

    exit_code = prepare_cache.main(
        [
            "--raw-interactions-path",
            str(raw_path),
            "--host-assembly-dir",
            str(host_dir),
            "--phage-fasta-dir",
            str(phage_dir),
            "--output-root",
            str(output_root),
            "--cache-dir",
            str(cache_dir),
            "--skip-host-assembly-resolution",
            "--holdout-fraction",
            "0.2",
            "--inner-val-fraction",
            "0.2",
            "--tl17-output-dir",
            str(tl17_output_dir),
        ]
    )
    assert exit_code == 0

    train_rows = read_csv_rows(cache_dir / "search_pairs" / prepare_cache.TRAIN_PAIR_TABLE_FILENAME)
    inner_val_rows = read_csv_rows(cache_dir / "search_pairs" / prepare_cache.INNER_VAL_PAIR_TABLE_FILENAME)
    assert train_rows
    assert inner_val_rows
    assert {row["split"] for row in train_rows} == {"train"}
    assert {row["split"] for row in inner_val_rows} == {"inner_val"}

    cache_manifest = json.loads((cache_dir / prepare_cache.CACHE_MANIFEST_FILENAME).read_text(encoding="utf-8"))
    assert set(cache_manifest["pair_tables"]) == {"train", "inner_val"}
    assert cache_manifest["feature_slots"]["host_defense"]["column_count"] == 1
    assert cache_manifest["feature_slots"]["phage_projection"]["column_count"] == 2
    assert cache_manifest["feature_slots"]["phage_stats"]["column_count"] == 4

    provenance = json.loads((cache_dir / prepare_cache.PROVENANCE_MANIFEST_FILENAME).read_text(encoding="utf-8"))
    assert provenance["build_mode"] == "raw_inputs_only"
    assert provenance["sealed_holdout"]["exported_to_search_cache"] is False
    assert provenance["sealed_holdout"]["retained_row_count"] > 0
    assert (
        provenance["search_workspace"]["feature_slots"]["host_surface"]["build_runtime"]["runtime_id"]
        == "test_fast_path"
    )
    assert provenance["search_workspace"]["feature_slots"]["host_surface"]["build_runtime"][
        "legacy_nhmmer_path_forbidden"
    ]
    assert provenance["search_workspace"]["feature_slots"]["phage_projection"]["build_runtime"]["projection_mode"] == (
        "batched_mmseqs_easy_search"
    )
    assert (
        provenance["search_workspace"]["feature_slots"]["phage_projection"]["reference_bank_provenance"]["task_id"]
        == "TL17"
    )

    host_index_rows = read_csv_rows(cache_dir / "feature_slots" / "host_surface" / prepare_cache.ENTITY_INDEX_FILENAME)
    exported_bacteria = {row["bacteria"] for row in host_index_rows}
    source_pair_rows = read_csv_rows(output_root / build_contract.PAIR_TABLE_FILENAME)
    holdout_bacteria = {row["bacteria"] for row in source_pair_rows if row["split"] == "holdout"}
    assert exported_bacteria
    # Feature slots include holdout entities (for AR09 replication embeddings);
    # only pair tables exclude holdout.
    assert holdout_bacteria.issubset(exported_bacteria)
    assert len(fast_path_calls) == 1
    assert fast_path_calls[0]["include_lps_core_type"] is False
    assert set(fast_path_calls[0]["assemblies"]) == exported_bacteria
    exported_phages = {
        row["phage"]
        for row in read_csv_rows(cache_dir / "feature_slots" / "phage_projection" / prepare_cache.ENTITY_INDEX_FILENAME)
    }
    assert len(phage_projection_calls) == 1
    assert set(phage_projection_calls[0]["phages"]) == exported_phages
    assert phage_projection_calls[0]["reference_fasta_path"] == tl17_output_dir / "tl17_rbp_reference_bank.faa"
    assert {call["phage"] for call in phage_stats_row_calls} == exported_phages

    feature_rows = read_csv_rows(cache_dir / "feature_slots" / "host_surface" / prepare_cache.SLOT_FEATURES_FILENAME)
    assert {key for key in feature_rows[0]} == {
        "bacteria",
        "host_surface__host_o_antigen_type",
        "host_surface__host_o_antigen_score",
        "host_surface__host_receptor_btub_score",
        "host_surface__host_capsule_profile_kpsc_score",
    }
    assert "host_surface__host_lps_core_type" not in feature_rows[0]

    slot_schema = json.loads(
        (cache_dir / "feature_slots" / "host_surface" / prepare_cache.SLOT_SCHEMA_FILENAME).read_text(encoding="utf-8")
    )
    assert slot_schema["reserved_feature_column_count"] == 4
    assert slot_schema["materialization"]["legacy_nhmmer_path_forbidden"] is True
    assert slot_schema["materialization"]["source_columns_dropped_for_autoresearch"] == ["host_lps_core_type"]

    schema = json.loads((cache_dir / prepare_cache.SCHEMA_MANIFEST_FILENAME).read_text(encoding="utf-8"))
    assert schema["feature_slots"]["host_defense"]["reserved_feature_columns"] == ["host_defense__AbiD"]
    assert schema["feature_slots"]["host_typing"]["reserved_feature_columns"] == [
        "host_typing__host_clermont_phylo",
        "host_typing__host_h_type",
        "host_typing__host_o_type",
        "host_typing__host_serotype",
        "host_typing__host_st_warwick",
    ]
    assert schema["feature_slots"]["host_stats"]["reserved_feature_columns"] == [
        "host_stats__host_gc_content",
        "host_stats__host_genome_length_nt",
        "host_stats__host_n50_contig_length_nt",
        "host_stats__host_sequence_record_count",
    ]
    assert schema["feature_slots"]["phage_projection"]["reserved_feature_columns"] == [
        "phage_projection__tl17_phage_rbp_family_11_percent_identity",
        "phage_projection__tl17_rbp_reference_hit_count",
    ]
    assert schema["feature_slots"]["phage_stats"]["reserved_feature_columns"] == [
        "phage_stats__phage_gc_content",
        "phage_stats__phage_genome_length_nt",
        "phage_stats__phage_n50_contig_length_nt",
        "phage_stats__phage_sequence_record_count",
    ]
    host_defense_rows = read_csv_rows(
        cache_dir / "feature_slots" / "host_defense" / prepare_cache.SLOT_FEATURES_FILENAME
    )
    assert {row["bacteria"] for row in host_defense_rows} == exported_bacteria

    host_typing_rows = read_csv_rows(cache_dir / "feature_slots" / "host_typing" / prepare_cache.SLOT_FEATURES_FILENAME)
    host_stats_rows = read_csv_rows(cache_dir / "feature_slots" / "host_stats" / prepare_cache.SLOT_FEATURES_FILENAME)
    phage_projection_rows = read_csv_rows(
        cache_dir / "feature_slots" / "phage_projection" / prepare_cache.SLOT_FEATURES_FILENAME
    )
    phage_stats_rows = read_csv_rows(cache_dir / "feature_slots" / "phage_stats" / prepare_cache.SLOT_FEATURES_FILENAME)
    typing_by_bacteria = {row["bacteria"]: row for row in host_typing_rows}
    stats_by_bacteria = {row["bacteria"]: row for row in host_stats_rows}
    projection_by_phage = {row["phage"]: row for row in phage_projection_rows}
    stats_by_phage = {row["phage"]: row for row in phage_stats_rows}
    joined_rows = []
    for row in train_rows:
        bacteria = row["bacteria"]
        phage = row["phage"]
        if row["retained_for_autoresearch"] != "1":
            continue
        joined_rows.append(
            {
                **row,
                **typing_by_bacteria[bacteria],
                **stats_by_bacteria[bacteria],
                **projection_by_phage[phage],
                **stats_by_phage[phage],
            }
        )
    assert joined_rows
    assert all("host_typing__host_serotype" in row for row in joined_rows)
    assert all("host_stats__host_gc_content" in row for row in joined_rows)
    assert all("phage_projection__tl17_rbp_reference_hit_count" in row for row in joined_rows)
    assert all("phage_stats__phage_gc_content" in row for row in joined_rows)

    host_typing_manifest = json.loads(
        (cache_dir / "feature_slots" / "host_typing" / prepare_cache.HOST_TYPING_BUILD_MANIFEST_FILENAME).read_text(
            encoding="utf-8"
        )
    )
    assert host_typing_manifest["guardrails"]["panel_metadata_used_for_feature_construction"] is False
    assert {entry["bacteria"] for entry in host_typing_manifest["runtime_caveats"]} == exported_bacteria

    phage_projection_manifest = json.loads(
        (
            cache_dir / "feature_slots" / "phage_projection" / prepare_cache.PHAGE_PROJECTION_BUILD_MANIFEST_FILENAME
        ).read_text(encoding="utf-8")
    )
    assert phage_projection_manifest["guardrails"]["batched_projection_path_used"] is True
    assert phage_projection_manifest["guardrails"]["checked_in_projection_csv_used"] is False
    assert phage_projection_manifest["reference_bank_provenance"]["reference_fasta_path"] == str(
        tl17_output_dir / "tl17_rbp_reference_bank.faa"
    )

    phage_stats_manifest = json.loads(
        (cache_dir / "feature_slots" / "phage_stats" / prepare_cache.PHAGE_STATS_BUILD_MANIFEST_FILENAME).read_text(
            encoding="utf-8"
        )
    )
    assert phage_stats_manifest["guardrails"]["per_phage_intermediate_artifacts_written"] is False

    phage_stats_schema = json.loads(
        (cache_dir / "feature_slots" / "phage_stats" / prepare_cache.SLOT_SCHEMA_FILENAME).read_text(encoding="utf-8")
    )
    assert phage_stats_schema["materialization"]["direct_feature_row_path_used"] is True


def test_validate_warm_cache_manifest_rejects_schema_mismatch(tmp_path: Path) -> None:
    schema = prepare_cache.build_top_level_schema_manifest()
    warm_cache_dir = tmp_path / "warm"
    warm_cache_dir.mkdir(parents=True)
    artifact_path = warm_cache_dir / "host_surface.csv"
    artifact_path.write_text("bacteria,host_surface__capsule_hits\nB1,1\n", encoding="utf-8")

    manifest_path = warm_cache_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "warm_cache_manifest_id": "fixture",
                "schema_manifest_id": "wrong_schema",
                "slot_artifacts": {
                    "host_surface": {
                        "path": "host_surface.csv",
                        "join_keys": ["bacteria"],
                        "column_family_prefix": "host_surface__",
                        "columns": ["host_surface__capsule_hits"],
                    }
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="schema mismatch"):
        prepare_cache.validate_warm_cache_manifest(manifest_path, schema_manifest=schema)


def test_validate_warm_cache_manifest_accepts_matching_descriptor(tmp_path: Path) -> None:
    schema = prepare_cache.build_top_level_schema_manifest()
    warm_cache_dir = tmp_path / "warm"
    warm_cache_dir.mkdir(parents=True)
    artifact_path = warm_cache_dir / "phage_stats.csv"
    artifact_path.write_text("phage,phage_stats__gc_content\nP1,0.5\n", encoding="utf-8")

    manifest_path = warm_cache_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "warm_cache_manifest_id": "fixture",
                "schema_manifest_id": prepare_cache.SCHEMA_MANIFEST_ID,
                "source_kind": "deploy_feature_csv_optional",
                "slot_artifacts": {
                    "phage_stats": {
                        "path": "phage_stats.csv",
                        "join_keys": ["phage"],
                        "column_family_prefix": "phage_stats__",
                        "columns": ["phage_stats__gc_content"],
                    }
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    validation = prepare_cache.validate_warm_cache_manifest(manifest_path, schema_manifest=schema)
    assert validation["schema_manifest_id"] == prepare_cache.SCHEMA_MANIFEST_ID
    assert validation["validated_slots"]["phage_stats"]["column_count"] == 1
    assert validation["validated_slots"]["phage_stats"]["columns"] == ["phage_stats__gc_content"]


def test_process_one_host_defense_cache_entry_forbids_worker_model_installs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    assembly_path = tmp_path / "B1.fna"
    assembly_path.write_text(">B1\nATGC\n", encoding="utf-8")
    calls: list[dict[str, object]] = []

    from lyzortx.pipeline.deployment_paired_features import derive_host_defense_features as host_defense_module

    def fake_derive_host_defense_features(
        assembly_path: Path,
        *,
        bacteria_id: str | None = None,
        output_dir: Path,
        panel_defense_subtypes_path: Path,
        models_dir: Path,
        workers: int,
        force_model_update: bool,
        model_install_mode: str,
        force_run: bool,
        preserve_raw: bool,
    ) -> dict[str, object]:
        calls.append(
            {
                "bacteria_id": bacteria_id,
                "output_dir": output_dir,
                "workers": workers,
                "force_model_update": force_model_update,
                "model_install_mode": model_install_mode,
            }
        )
        return {"feature_row": {"bacteria": bacteria_id or assembly_path.stem}}

    monkeypatch.setattr(host_defense_module, "derive_host_defense_features", fake_derive_host_defense_features)

    bacteria_id, ok, message = prepare_cache._process_one_host_defense_cache_entry(
        assembly_path,
        "B1",
        tmp_path / "per_host",
        Path("data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv"),
        tmp_path / "models",
        False,
        False,
    )

    assert (bacteria_id, ok, message) == ("B1", True, "ok")
    assert calls == [
        {
            "bacteria_id": "B1",
            "output_dir": tmp_path / "per_host" / "B1",
            "workers": 1,
            "force_model_update": False,
            "model_install_mode": "forbid",
        }
    ]


def test_build_host_defense_slot_artifact_supports_reaggregation_without_worker_rerun(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cache_dir = tmp_path / "cache"
    output_root = tmp_path / "outputs"
    per_host_output_dir = output_root / "host_defense"
    for bacteria in ("B1", "B2"):
        host_dir = per_host_output_dir / bacteria
        host_dir.mkdir(parents=True, exist_ok=True)
        (host_dir / prepare_cache.PER_HOST_COUNTS_FILENAME).write_text(
            f"bacteria,AbiD,RM_Type_I\n{bacteria},1,0\n",
            encoding="utf-8",
        )

    def fail_if_called(*args: object, **kwargs: object) -> str:
        raise AssertionError("ensure_defense_finder_models should not run in aggregate-only mode")

    monkeypatch.setattr(prepare_cache, "ensure_defense_finder_models", fail_if_called)

    args = prepare_cache.parse_args(
        [
            "--output-root",
            str(output_root),
            "--cache-dir",
            str(cache_dir),
            "--skip-host-assembly-resolution",
            "--host-defense-output-dir",
            str(per_host_output_dir),
            "--host-defense-aggregate-only",
        ]
    )
    split_rows = {
        "train": [
            {
                "bacteria": "B1",
                "host_fasta_path": str(tmp_path / "hosts" / "B1.fna"),
                "retained_for_autoresearch": "1",
            }
        ],
        "inner_val": [
            {
                "bacteria": "B2",
                "host_fasta_path": str(tmp_path / "hosts" / "B2.fna"),
                "retained_for_autoresearch": "1",
            }
        ],
    }

    summary = prepare_cache._build_host_defense_slot_artifact(
        args=args,
        cache_dir=cache_dir,
        split_rows=split_rows,
    )

    assert summary["column_count"] > 0
    artifact_rows = read_csv_rows(cache_dir / "feature_slots" / "host_defense" / prepare_cache.SLOT_FEATURES_FILENAME)
    assert [row["bacteria"] for row in artifact_rows] == ["B1", "B2"]
    assert all(row["host_defense__AbiD"] == "1" for row in artifact_rows)
    assert all(row["host_defense__RM_Type_I"] == "0" for row in artifact_rows)
    assert all(column == "bacteria" or column.startswith("host_defense__") for column in artifact_rows[0])


def test_build_host_defense_slot_artifact_wraps_process_pool_failures(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cache_dir = tmp_path / "cache"
    output_root = tmp_path / "outputs"
    process_error = OSError("spawn failed")
    logged_errors: list[str] = []

    class FakeFuture:
        def result(self) -> tuple[str, bool, str]:
            raise process_error

    class FakeProcessPoolExecutor:
        def __init__(self, *, max_workers: int) -> None:
            self.max_workers = max_workers

        def __enter__(self) -> FakeProcessPoolExecutor:
            return self

        def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
            return None

        def submit(self, *args: object, **kwargs: object) -> FakeFuture:
            return FakeFuture()

    def fake_as_completed(futures: dict[FakeFuture, str]) -> list[FakeFuture]:
        return list(futures)

    def fake_ensure_defense_finder_models(*args: object, **kwargs: object) -> str:
        return "already_present"

    def fail_if_called(*args: object, **kwargs: object) -> None:
        raise AssertionError("aggregate_host_defense_csvs should not run when a worker process fails")

    monkeypatch.setattr(prepare_cache, "ProcessPoolExecutor", FakeProcessPoolExecutor)
    monkeypatch.setattr(prepare_cache, "as_completed", fake_as_completed)
    monkeypatch.setattr(prepare_cache, "ensure_defense_finder_models", fake_ensure_defense_finder_models)
    monkeypatch.setattr(prepare_cache, "aggregate_host_defense_csvs", fail_if_called)
    monkeypatch.setattr(
        prepare_cache.LOGGER,
        "error",
        lambda message, *args: logged_errors.append(message % args),
    )

    args = prepare_cache.parse_args(
        [
            "--output-root",
            str(output_root),
            "--cache-dir",
            str(cache_dir),
            "--skip-host-assembly-resolution",
        ]
    )
    split_rows = {
        "train": [
            {
                "bacteria": "B1",
                "host_fasta_path": str(tmp_path / "hosts" / "B1.fna"),
                "retained_for_autoresearch": "1",
            }
        ],
        "inner_val": [],
    }

    with pytest.raises(RuntimeError, match="AR03 host-defense cache build failed for B1"):
        prepare_cache._build_host_defense_slot_artifact(
            args=args,
            cache_dir=cache_dir,
            split_rows=split_rows,
        )

    assert logged_errors == ["AR03 host defense failed for B1: worker process error: spawn failed"]
