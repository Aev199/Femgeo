#!/usr/bin/env bash
# Изменения OpenCS в перенесённых файлах между закреплённым коммитом и новым.
#   scripts/upstream_diff.sh [новый SHA или ветка, по умолчанию master] [--patch]
# Клон OpenCS хранится вне репозитория: ${FEMGEO_UPSTREAM_DIR:-$HOME/.cache/femgeo/OpenCS}.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
lock="$root/upstream/opencs.lock"
target="${1:-master}"
patch="${2:-}"
repo=$(sed -n 's/^repository = //p' "$lock")
old=$(sed -n 's/^commit = //p' "$lock")
dir="${FEMGEO_UPSTREAM_DIR:-$HOME/.cache/femgeo/OpenCS}"
if [ ! -d "$dir/.git" ]; then
  mkdir -p "$(dirname "$dir")"
  git clone --quiet "$repo" "$dir"
fi
git -C "$dir" fetch --quiet origin
new=$(git -C "$dir" rev-parse --verify "origin/$target" 2>/dev/null || git -C "$dir" rev-parse --verify "$target")
mapfile -t files < <(grep -v '^#' "$lock" | grep ' = ' | grep -v '^repository\|^commit' \
  | sed 's/^[^=]*= //' | tr ',' '\n' | sed 's/^ *//; s/ *$//' | sort -u)
echo "OpenCS: $old -> $new"
changed=$(git -C "$dir" diff --name-only "$old" "$new" -- "${files[@]}" || true)
if [ -z "$changed" ]; then
  echo "Перенесённые файлы не менялись."
  exit 0
fi
echo "Изменены перенесённые файлы:"
echo "$changed" | sed 's/^/  /'
echo
echo "Наши файлы, которые нужно сверить:"
while IFS= read -r f; do
  grep -v '^#' "$lock" | grep ' = ' | grep -F "$f" | sed 's/ = .*//; s/^/  /'
done <<< "$changed" | sort -u
if [ "$patch" = "--patch" ]; then
  git -C "$dir" diff "$old" "$new" -- $changed
fi
