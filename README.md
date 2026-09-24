# Rare Variant Explorer

A compact portfolio project for exploring **public, simple GRCh38 variants**. It connects a VCF, Ensembl VEP, PostgreSQL, an optional Human Phenotype Ontology (HPO) import, a FastAPI service, and a small JavaScript browser UI. It is a **research demo**, not a clinical interpretation or patient matching system.

```text
public GRCh38 VCF ──> Ensembl VEP REST or local VEP Docker ──> JSONL
                                                              │
HPO genes_to_phenotype.txt ──> gene / disease / HPO tables ──> PostgreSQL
                                                              │
                                               FastAPI ──> browser UI
```

## Quick start (small VCF, live Ensembl REST)

Requires Docker Compose and internet access to Ensembl REST. These example variants are from the [official Ensembl VEP region POST documentation](https://rest.ensembl.org/documentation/info/vep_region_post); the VCF declares GRCh38. **Do not substitute patient VCFs for the REST demonstration.**

```bash
docker compose up -d --build
docker compose run --rm app python -m explorer.cli annotate-rest data/example.vcf --output data/vep.jsonl
docker compose run --rm app python -m explorer.cli import-vep data/example.vcf data/vep.jsonl --source ensembl-rest --version 'Ensembl REST release 116'
```

Use the **actual release number printed** by `annotate-rest` in the third command; 116 above is an example. Open <http://localhost:8001>. API documentation: <http://localhost:8001/docs>. `GET /api/variants?gene=BRCA2` searches symbols and `GET /api/variants/{id}` shows a complete record. The Docker host port is 8001 to avoid conflicts with other local FastAPI services. If Ensembl is unavailable or rate limited, run the local Docker VEP path below instead. No invented annotation results are bundled.

### More useful public example: CFTR p.Gly551Asp (G551D)

The separate `data/cftr_g551d.vcf` contains a public GRCh38 single nucleotide allele corresponding to [ClinVar CFTR p.Gly551Asp](https://www.ncbi.nlm.nih.gov/clinvar/variation/7120/). Its VCF coordinates and bases are `7:117587806 G>A`. The condition-specific [ClinVar RCV record](https://www.ncbi.nlm.nih.gov/clinvar/RCV000007540/) is bundled as a clearly dated snapshot in `data/clinvar_g551d.json`; it is imported **separately** and does not imply that VEP produces pathogenicity assertions. The older p.Phe508del indel sample remains in `data/cftr_f508del.vcf`, but the public Ensembl REST endpoint returned HTTP 500 for it on 2026-09-24; use G551D for this walkthrough.

After the first quick start, run:

```bash
docker compose run --rm app python -m explorer.cli annotate-rest data/cftr_g551d.vcf --output data/cftr_g551d_vep.jsonl
docker compose run --rm app python -m explorer.cli import-vep data/cftr_g551d.vcf data/cftr_g551d_vep.jsonl --source ensembl-rest --version 'Ensembl REST release 116'
docker compose run --rm app python -m explorer.cli import-clinvar-example data/clinvar_g551d.json
curl -L --fail 'https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt' -o data/genes_to_phenotype.txt
docker compose run --rm app python -m explorer.cli import-hpo data/genes_to_phenotype.txt
```

Again, use the release actually printed by `annotate-rest`. Search **CFTR** in the browser and click its row. ClinVar's `Pathogenic` classification is a source assertion for **cystic fibrosis (OMIM:219700)**, with a record accession and review status. The HPO panel shows CFTR–phenotype associations **restricted to the same disease identifier**, with the original HPO source fingerprint. A gene–disease HPO association is still not evidence that an individual carrying this allele has that feature. The live gnomAD genomes value may still be missing from Ensembl REST; no value is imputed.

### macOS Ventura 13 / Intel: run without Docker

New Docker Desktop versions require macOS 14+. On Ventura, install [Postgres.app with PostgreSQL 16](https://postgresapp.com/downloads.html) (the universal macOS build), move it to Applications, launch it and click **Initialize**. Use the Python 3.14 you already have, or any installed Python 3.12+. From Terminal, in the extracted project folder:

```bash
/Applications/Postgres.app/Contents/Versions/latest/bin/createdb explorer
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export DATABASE_URL="postgresql://$(whoami)@localhost:5432/explorer"
python -m explorer.cli init
python -m explorer.cli annotate-rest data/example.vcf --output data/vep.jsonl
python -m explorer.cli import-vep data/example.vcf data/vep.jsonl --source ensembl-rest --version 'Ensembl REST release 116'
python -m uvicorn explorer.api:app --host 127.0.0.1 --port 8000
```

Replace `116` in the import command with the actual release printed by the annotate command. Open <http://127.0.0.1:8000>. Leave that Terminal open while using the app; press **Control-C** to stop the web server. Postgres.app must remain running. Its default local server uses port 5432 and creates a role matching your macOS user; change the `DATABASE_URL` if you configured a different port, role or password. An already existing `explorer` database makes `createdb` report an error; you can proceed to `init` in that case.

### HPO associations (optional)

Download the [official `genes_to_phenotype.txt`](https://obophenotype.github.io/human-phenotype-ontology/annotations/genes_to_phenotype/) file, then import it after VEP:

```bash
curl -L --fail 'https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt' -o data/genes_to_phenotype.txt
docker compose run --rm app python -m explorer.cli import-hpo data/genes_to_phenotype.txt
```

The importer reads only associations for genes found in imported variants and records a SHA-256 of the HPO download. A gene–HPO–disease row is **not** a variant–phenotype observation; symbol matching is a demo simplification. An empty panel can mean no transcript gene symbol, no matching HPO entry, or no HPO import. Imported HPO terms are attributed to the [Human Phenotype Ontology](https://obophenotype.github.io/human-phenotype-ontology/); the source SHA-256 identifies the downloaded release bytes.

### Local VEP Docker path

For reproducible offline annotation, install the **Ensembl release 116 GRCh38 VEP cache** in a local directory and set `VEP_CACHE` to its absolute path. Consult the [VEP download and cache instructions](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) for the matching cache; it is large and is intentionally excluded from this repository. The script uses `ensemblorg/ensembl-vep:release_116`, `--offline`, `--af_gnomadg`, JSON output, and the same example VCF.

```bash
export VEP_CACHE=/absolute/path/to/vep_cache
bash scripts/annotate_docker.sh
docker compose run --rm app python -m explorer.cli import-vep data/example.vcf data/vep.jsonl --source vep-docker --version 'Ensembl VEP release 116; GRCh38 cache 116'
```

An absent gnomAD genome frequency stays **null**. A `0` frequency is kept as zero. VEP REST may not expose the same cache population fields as local VEP; this app never fills missing frequencies with an estimate. Population frequency is not a pathogenicity classification.

## Design decisions and constraints

- `schema.sql` separates annotation runs, alleles, ClinVar condition assertions, HPO terms and gene–disease–phenotype associations. It stores the full VEP JSON to avoid discarding alternative transcript consequences. The UI's selected transcript prefers MANE Select, then canonical, but the variant's most severe consequence is VEP's site-level field, which can describe a **different transcript**.
- Importing the same VCF, source and source-version twice updates the same run. The SHA-256 of VCF bytes is stored. Import fails atomically if VEP returns an extra, duplicate or missing allele; it checks assembly when the annotation specifies one.
- `source_version` is supplied by the operator. Check it against the REST release reported during annotation or the installed Docker image/cache; the importer cannot independently authenticate the declared version. The REST output is saved to JSONL for inspection and repeat import.
- The bundled ClinVar records are **2026-09-24 snapshots**, not a synced feed. Check the linked live condition record before using a classification in any other context. Import refuses to attach an assertion to an allele unless that exact GRCh38 allele and gene have already been annotated and stored.
- Input is capped at 1,000 simple A/C/G/T alleles, with required `##reference=GRCh38`. It does not validate reference bases against a FASTA, left normalize indels, handle structural variants, interpret genotypes, classify variants, handle patient consent, or verify gene symbol aliases. Avoid using the application with protected clinical data.
- PostgreSQL is mapped only to local port `5433`; the Docker web interface is mapped only to local `8001` (the native Mac setup uses `8000`). The Compose password is for a local demo only. Add authentication, TLS, credential management, migrations and audit controls before any shared deployment.
- This is deliberately plain JavaScript rather than a framework: query, selection, pagination and detail rendering are visible without a build step. It demonstrates browser-side development, while a production Vue interface remains future work.

## Verification

The parser and annotation tests use the Python standard library; they need no database or internet:

```bash
python3 -m unittest discover -s tests -v
node --check static/app.js
docker compose config
```

With Docker available, run the quick start and check `/api/health`, the displayed variant rows, full JSON, and (after an HPO import) the gene association panel. The sample VCF and HPO mappings contain **public data only**. For an independent variant biology check, inspect the recorded VEP JSON and the source databases; this application is not validated for clinical decision making.

## Sources

- [Ensembl VEP REST region POST format](https://rest.ensembl.org/documentation/info/vep_region_post) and [VEP download / Docker instructions](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)
- [Ensembl VEP CLI options, including `--af_gnomadg`](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)
- [HPO gene-to-phenotype file format](https://obophenotype.github.io/human-phenotype-ontology/annotations/genes_to_phenotype/)
- [ClinVar variation record VCV000007120](https://www.ncbi.nlm.nih.gov/clinvar/variation/7120/) and [condition record RCV000007540](https://www.ncbi.nlm.nih.gov/clinvar/RCV000007540/)
