#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VEP_CACHE:?Set VEP_CACHE to the directory holding the human GRCh38 release 116 cache}"
if [[ ! -f data/example.vcf ]]; then
  echo 'Missing data/example.vcf' >&2
  exit 1
fi
docker run --rm \
  -v "$PWD/data:/work" \
  -v "$VEP_CACHE:/cache:ro" \
  ensemblorg/ensembl-vep:release_116 \
  vep --input_file /work/example.vcf --output_file /work/vep.jsonl \
  --format vcf --json --no_stats --force_overwrite \
  --assembly GRCh38 --cache --offline --dir_cache /cache \
  --symbol --canonical --mane --af_gnomadg
echo 'VEP annotations saved to data/vep.jsonl'
