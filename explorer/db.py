import os
from pathlib import Path

import psycopg
from psycopg.rows import dict_row
from psycopg.types.json import Jsonb

from explorer.processing import annotations_for_records, load_clinvar_example, parse_hpo_gene_associations, parse_vcf, sha256


def connect():
    return psycopg.connect(os.environ.get("DATABASE_URL", "postgresql://explorer:localdemo@localhost:5433/explorer"), row_factory=dict_row)


def init():
    with connect() as conn:
        conn.execute(Path(__file__).resolve().parent.parent.joinpath("schema.sql").read_text())


def import_annotations(vcf, objects, source, version):
    records = parse_vcf(vcf)
    parsed = annotations_for_records(records, objects)
    with connect() as conn:
        run_id = conn.execute("""
            INSERT INTO annotation_runs (input_sha256,source,source_version,assembly)
            VALUES (%s,%s,%s,'GRCh38')
            ON CONFLICT (input_sha256,source,source_version)
            DO UPDATE SET source=EXCLUDED.source RETURNING id
        """, (sha256(vcf), source, version)).fetchone()["id"]
        for (chrom, pos, ref, alt), summary, raw in parsed:
            conn.execute("""
                INSERT INTO variants (run_id,chrom,pos,ref,alt,gene_id,gene_symbol,
                  transcript_id,consequence,gnomadg_af,vep_json)
                VALUES (%(run_id)s,%(chrom)s,%(pos)s,%(ref)s,%(alt)s,%(gene_id)s,
                  %(gene_symbol)s,%(transcript_id)s,%(consequence)s,%(gnomadg_af)s,%(raw)s)
                ON CONFLICT (run_id,chrom,pos,ref,alt) DO UPDATE SET
                  gene_id=EXCLUDED.gene_id, gene_symbol=EXCLUDED.gene_symbol,
                  transcript_id=EXCLUDED.transcript_id,
                  consequence=EXCLUDED.consequence, gnomadg_af=EXCLUDED.gnomadg_af,
                  vep_json=EXCLUDED.vep_json
            """, {"run_id": run_id, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                  **summary, "raw": Jsonb(raw)})
    return run_id, len(parsed)


def import_hpo(path):
    with connect() as conn:
        genes = {r["gene_symbol"] for r in conn.execute("SELECT DISTINCT gene_symbol FROM variants WHERE gene_symbol IS NOT NULL")}
    associations = parse_hpo_gene_associations(path, gene_filter=genes)
    digest = sha256(path)
    with connect() as conn:
        for gene, hpo_id, name, disease in associations:
            conn.execute("""INSERT INTO hpo_terms (hpo_id,name) VALUES (%s,%s)
                            ON CONFLICT (hpo_id) DO UPDATE SET name=EXCLUDED.name""", (hpo_id, name))
            conn.execute("""INSERT INTO gene_phenotypes (gene_symbol,hpo_id,disease_id,source_sha256)
                            VALUES (%s,%s,%s,%s) ON CONFLICT DO NOTHING""", (gene, hpo_id, disease, digest))
    return len(associations)


def import_clinvar_example(path):
    assertions = load_clinvar_example(path)
    with connect() as conn:
        for row in assertions:
            match = conn.execute("""SELECT 1 FROM variants WHERE chrom=%s AND pos=%s
                AND ref=%s AND alt=%s AND gene_symbol=%s LIMIT 1""",
                (row["chrom"], row["pos"], row["ref"], row["alt"], row["gene_symbol"])).fetchone()
            if not match:
                raise ValueError(f"Annotate and import {row['gene_symbol']} {row['chrom']}:{row['pos']} first")
            conn.execute("""INSERT INTO clinvar_assertions
                (assembly,chrom,pos,ref,alt,gene_symbol,condition,disease_id,
                 classification,review_status,record_accession,source_url,snapshot_date)
                VALUES (%(assembly)s,%(chrom)s,%(pos)s,%(ref)s,%(alt)s,%(gene_symbol)s,
                        %(condition)s,%(disease_id)s,%(classification)s,%(review_status)s,
                        %(record_accession)s,%(source_url)s,%(snapshot_date)s)
                ON CONFLICT (assembly,chrom,pos,ref,alt,record_accession) DO UPDATE SET
                    classification=EXCLUDED.classification, review_status=EXCLUDED.review_status,
                    snapshot_date=EXCLUDED.snapshot_date, source_url=EXCLUDED.source_url""", row)
    return len(assertions)
