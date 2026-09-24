CREATE TABLE IF NOT EXISTS annotation_runs (
    id BIGSERIAL PRIMARY KEY,
    input_sha256 CHAR(64) NOT NULL,
    source TEXT NOT NULL,
    source_version TEXT NOT NULL,
    assembly TEXT NOT NULL CHECK (assembly = 'GRCh38'),
    annotated_at TIMESTAMPTZ NOT NULL DEFAULT now(),
    UNIQUE (input_sha256, source, source_version)
);

CREATE TABLE IF NOT EXISTS variants (
    id BIGSERIAL PRIMARY KEY,
    run_id BIGINT NOT NULL REFERENCES annotation_runs(id) ON DELETE CASCADE,
    chrom TEXT NOT NULL,
    pos BIGINT NOT NULL CHECK (pos > 0),
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    gene_id TEXT,
    gene_symbol TEXT,
    transcript_id TEXT,
    consequence TEXT,
    gnomadg_af DOUBLE PRECISION CHECK (gnomadg_af BETWEEN 0 AND 1),
    vep_json JSONB NOT NULL,
    UNIQUE (run_id, chrom, pos, ref, alt)
);
CREATE INDEX IF NOT EXISTS variants_gene_idx ON variants (gene_symbol);
CREATE INDEX IF NOT EXISTS variants_region_idx ON variants (chrom, pos);

CREATE TABLE IF NOT EXISTS hpo_terms (
    hpo_id TEXT PRIMARY KEY,
    name TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS gene_phenotypes (
    gene_symbol TEXT NOT NULL,
    hpo_id TEXT NOT NULL REFERENCES hpo_terms(hpo_id),
    disease_id TEXT NOT NULL,
    source_sha256 CHAR(64) NOT NULL,
    PRIMARY KEY (gene_symbol, hpo_id, disease_id, source_sha256)
);
CREATE INDEX IF NOT EXISTS gene_phenotypes_gene_idx ON gene_phenotypes (gene_symbol);

CREATE TABLE IF NOT EXISTS clinvar_assertions (
    id BIGSERIAL PRIMARY KEY,
    assembly TEXT NOT NULL CHECK (assembly='GRCh38'),
    chrom TEXT NOT NULL,
    pos BIGINT NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    gene_symbol TEXT NOT NULL,
    condition TEXT NOT NULL,
    disease_id TEXT NOT NULL,
    classification TEXT NOT NULL,
    review_status TEXT NOT NULL,
    record_accession TEXT NOT NULL,
    source_url TEXT NOT NULL,
    snapshot_date DATE NOT NULL,
    UNIQUE (assembly,chrom,pos,ref,alt,record_accession)
);
CREATE INDEX IF NOT EXISTS clinvar_assertions_allele_idx
    ON clinvar_assertions (chrom,pos,ref,alt);
