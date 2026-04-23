---
name: worktree
description: >-
  Git worktree workflow for Claude Code and Codex — creating, using, and cleaning up
  isolated worktrees for parallel development.
user-invocable: true
allowed-tools: Bash, Read, Glob
argument-hint: "[create <name> | list | cleanup [--dry-run]]"
---

# Git Worktree Workflow

Git worktrees let you check out multiple branches simultaneously in separate directories,
sharing the same repository history. They are ideal for:

- Working on a feature branch without disturbing the main checkout
- Running parallel agent sessions on different tasks
- Isolating experimental changes that may be discarded

## Canonical docs

- **Claude Code**: https://code.claude.com/docs/en/common-workflows#run-parallel-claude-code-sessions-with-git-worktrees
- **Codex**: https://developers.openai.com/codex/app/worktrees

## Worktree storage locations

| Tool        | Default location                                              |
|-------------|---------------------------------------------------------------|
| Claude Code | `<repo>/.claude/worktrees/<name>/`                            |
| Codex       | `$CODEX_HOME/worktrees/` (defaults to `~/.codex/worktrees/`) |

Both are gitignored and should stay that way.

## Subcommands

Route on the first word of `$ARGUMENTS`:

### `create <name>`

1. Run `git fetch origin` to get the latest remote state.
2. Create a worktree branching from `origin/main`:
   ```bash
   git worktree add .claude/worktrees/<name> -b worktree-<name> origin/main
   ```
3. Report the path and branch name so the user can `cd` into it or launch a session with
   `claude --worktree <name>` / `codex --worktree <name>`.

### `list`

1. Run `git worktree list` and display the output.
2. Also check for directories in `.claude/worktrees/` and `$CODEX_HOME/worktrees/`
   (or `~/.codex/worktrees/`) that may not be registered (orphaned directories).

### `cleanup [--dry-run]`

1. Run `git fetch --prune` to update remote tracking info.
2. Run `git worktree list --porcelain` to enumerate all worktrees.
3. For each worktree (excluding the main working tree):
   a. Identify its branch.
   b. Check if the branch still exists on the remote (`git ls-remote --heads origin <branch>`).
   c. Check if the worktree has uncommitted changes (`git -C <path> status --porcelain`).
   d. Classify as: **safe to remove** (remote branch gone, no local changes),
      **has local changes** (warn user), or **still active** (remote branch exists, skip).
4. Present a summary table showing each worktree, its branch, status, and recommended action.
5. If `--dry-run` is present, stop after the summary.
6. Otherwise, ask the user for confirmation, then:
   - Remove safe worktrees with `git worktree remove <path>`
   - Delete their local branches with `git branch -d <branch>`
   - Run `git worktree prune`
7. Report what was cleaned up.

### No arguments / help

Print a short usage summary of the subcommands above.

## Safety rules

- Never remove a worktree that has uncommitted changes without explicit user confirmation.
- Never remove the main working tree.
- Always show what will be done before doing it.
