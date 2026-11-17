#!/usr/bin/env bash
# Simple fake blastn wrapper for testing. Writes its argv to a known file and
# creates the -out file if provided.

INVOCATION_LOG="/tmp/fake_blastn_invocation.txt"
echo "$@" > "$INVOCATION_LOG"

OUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -out)
      OUT="$2"
      shift 2
      ;;
    -out=*)
      OUT="${1#*=}"
      shift
      ;;
    *)
      shift
      ;;
  esac
done

if [[ -n "$OUT" ]]; then
  mkdir -p "$(dirname "$OUT")" 2>/dev/null || true
  # write a minimal CSV header so pandas.read_csv can parse an empty-like file
  echo "" > "$OUT"
fi

exit 0
