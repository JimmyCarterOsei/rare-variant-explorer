from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from explorer.db import connect, init

ROOT = Path(__file__).resolve().parent.parent
app = FastAPI(title="Rare Variant Explorer", description="Public GRCh38 variants; research demo, no clinical classification")
app.mount("/static", StaticFiles(directory=ROOT / "static"), name="static")


@app.on_event("startup")
def startup():
    init()


@app.get("/", include_in_schema=False)
def home():
    return FileResponse(ROOT / "static/index.html")


@app.get("/api/health")
def health():
    with connect() as conn:
        conn.execute("SELECT 1")
    return {"status": "ok"}


@app.get("/api/variants")
def variants(gene: str | None = Query(default=None, max_length=40),
             limit: int = Query(default=50, ge=1, le=100),
             offset: int = Query(default=0, ge=0)):
    with connect() as conn:
        rows = conn.execute("""
            SELECT v.id,v.chrom,v.pos,v.ref,v.alt,v.gene_symbol,v.consequence,
                   v.gnomadg_af,r.source,r.source_version,r.assembly,r.annotated_at
            FROM variants v JOIN annotation_runs r ON r.id=v.run_id
            WHERE (%s::text IS NULL OR v.gene_symbol ILIKE %s)
            ORDER BY v.id DESC LIMIT %s OFFSET %s
        """, (gene, f"%{gene or ''}%", limit, offset)).fetchall()
    return {"items": rows, "limit": limit, "offset": offset}


@app.get("/api/variants/{variant_id}")
def variant(variant_id: int):
    with connect() as conn:
        record = conn.execute("""
            SELECT v.*,r.source,r.source_version,r.assembly,r.input_sha256,r.annotated_at
            FROM variants v JOIN annotation_runs r ON r.id=v.run_id WHERE v.id=%s
        """, (variant_id,)).fetchone()
        if record is None:
            raise HTTPException(404, "Variant not found")
        assertions = conn.execute("""
            SELECT condition,disease_id,classification,review_status,
                   record_accession,source_url,snapshot_date
            FROM clinvar_assertions
            WHERE assembly='GRCh38' AND chrom=%s AND pos=%s AND ref=%s AND alt=%s
            ORDER BY disease_id,record_accession
        """, (record["chrom"],record["pos"],record["ref"],record["alt"])).fetchall()
        conditions = sorted({a["disease_id"] for a in assertions})
        if record["gene_symbol"]:
            sql = """SELECT DISTINCT p.hpo_id,t.name,p.disease_id,p.source_sha256
                FROM gene_phenotypes p JOIN hpo_terms t ON t.hpo_id=p.hpo_id
                WHERE p.gene_symbol=%s"""
            params = [record["gene_symbol"]]
            if conditions:
                sql += " AND p.disease_id=ANY(%s::text[])"
                params.append(conditions)
            sql += " ORDER BY p.disease_id,p.hpo_id LIMIT 250"
            associations = conn.execute(sql, params).fetchall()
        else:
            associations = []
    record["clinvar_assertions"] = assertions
    record["hpo_condition_ids"] = conditions
    record["gene_phenotypes"] = associations
    return record
