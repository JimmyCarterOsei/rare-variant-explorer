"""Explicit, reproducible import commands (read-only annotation over public VCF)."""
import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

from explorer.processing import annotations_for_records, load_jsonl, parse_vcf


def fetch_rest(vcf, output):
    import httpx

    records = parse_vcf(vcf)
    url = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    objects = []
    with httpx.Client(timeout=90, headers=headers) as client:
        release = client.get(url + "/info/data").raise_for_status().json()
        for i in range(0, len(records), 100):
            chunk = records[i:i+100]
            for attempt in range(3):
                response = client.post(url + "/vep/human/region",
                                       params={"symbol": 1, "canonical": 1, "mane": 1},
                                       json={"variants": [row[4] for row in chunk]})
                if response.status_code not in (500, 502, 503, 504) or attempt == 2:
                    break
                time.sleep(2 ** attempt)
            if response.status_code >= 400:
                raise RuntimeError(f"Ensembl VEP HTTP {response.status_code}: {response.text[:400]}")
            response.raise_for_status()
            objects.extend(response.json())
    annotations_for_records(records, objects)  # Validate before saving.
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("".join(json.dumps(o, sort_keys=True) + "\n" for o in objects))
    version = "Ensembl REST release " + ",".join(str(x) for x in release.get("releases", []))
    return version, datetime.now(timezone.utc).isoformat()


def main():
    parser = argparse.ArgumentParser(description="Public variant exploration, GRCh38 only")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("init")
    rest = commands.add_parser("annotate-rest")
    rest.add_argument("vcf")
    rest.add_argument("--output", default="data/vep.jsonl")
    imported = commands.add_parser("import-vep")
    imported.add_argument("vcf")
    imported.add_argument("jsonl")
    imported.add_argument("--source", required=True, choices=["ensembl-rest", "vep-docker"])
    imported.add_argument("--version", required=True, help="For example: Ensembl 116, VEP release_116")
    hpo = commands.add_parser("import-hpo")
    hpo.add_argument("associations", help="Official genes_to_phenotype.txt")
    clinvar = commands.add_parser("import-clinvar-example")
    clinvar.add_argument("snapshot", help="Verified public condition assertion JSON")
    args = parser.parse_args()
    if args.command == "init":
        from explorer.db import init
        init()
        print("Database schema ready")
    elif args.command == "annotate-rest":
        version, timestamp = fetch_rest(args.vcf, args.output)
        print(f"Saved VEP JSON to {args.output}; {version}; fetched {timestamp}")
        print("Run import-vep with --source ensembl-rest --version and the release shown above")
    elif args.command == "import-vep":
        from explorer.db import init, import_annotations
        init()
        run_id, count = import_annotations(args.vcf, load_jsonl(args.jsonl), args.source, args.version)
        print(f"Imported {count} alleles to run {run_id}")
    elif args.command == "import-hpo":
        from explorer.db import init, import_hpo
        init()
        print(f"Imported {import_hpo(args.associations)} gene–HPO–disease associations")
    elif args.command == "import-clinvar-example":
        from explorer.db import import_clinvar_example, init
        init()
        print(f"Imported {import_clinvar_example(args.snapshot)} ClinVar condition assertions")


if __name__ == "__main__":
    main()
