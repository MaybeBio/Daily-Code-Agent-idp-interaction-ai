#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILLS_DIR="${HOME}/.claude/skills"

linked=0
skipped=0

for skill_dir in "$REPO_DIR"/*/; do
    skill_name="$(basename "$skill_dir")"
    [[ -f "${skill_dir}SKILL.md" ]] || continue

    target="${SKILLS_DIR}/${skill_name}"

    # Collision: exists and is NOT already a symlink to this repo
    if [[ -e "$target" && ! -L "$target" ]]; then
        echo "SKIP (conflict): $skill_name — real directory already exists at $target"
        skipped=$((skipped + 1))
        continue
    fi

    ln -sfn "$skill_dir" "$target"
    echo "Linked: $skill_name"
    linked=$((linked + 1))
done

echo "Done. $linked linked, $skipped skipped."
