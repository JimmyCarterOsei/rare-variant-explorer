"""Small, strict input and annotation adapters. No clinical classification."""
import hashlib
import json
import re
from datetime import date
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_vcf(path, limit=1000):
    """Accept simple GRCh38 VCF sites; split multiallelic sites for VEP.

    Genotypes are deliberately ignored. This is a public variant browser.
    """
    records = []
    saw_header = False
    assembly = False
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##reference=GRCh38"):
                assembly = True
            if line.startswith("#CHROM\t"):
                saw_header = True
                continue
            if line.startswith("#"):
                continue
            if not saw_header or not assembly:
                raise ValueError("VCF needs #CHROM header and ##reference=GRCh38")
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                raise ValueError("VCF record must contain at least 8 columns")
            chrom, position, identifier, ref, alts = parts[:5]
            chrom = chrom.removeprefix("chr")
            if chrom not in [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]:
                raise ValueError(f"Unsupported chromosome: {chrom}")
            pos = int(position)
            if pos < 1 or not ref or ref == ".":
                raise ValueError("Invalid position/reference")
            for alt in alts.split(","):
                if any(c not in "ACGT" for c in (ref + alt).upper()) or not alt:
                    raise ValueError("Only simple A/C/G/T alleles are supported")
                # VEP identifies its input with the original eight VCF columns.
                records.append((chrom, pos, ref.upper(), alt.upper(),
                                f"{chrom}\t{pos}\t{identifier}\t{ref}\t{alt}\t" + "\t".join(parts[5:8])))
                if len(records) > limit:
                    raise ValueError(f"Limit is {limit} alleles")
    if not records:
        raise ValueError("VCF has no variants")
    return records


def annotations_for_records(records, vep_objects):
    """Match by VEP input, fail on missing/extra/ambiguous responses."""
    expected = {(c, p, r, a): raw for c, p, r, a, raw in records}
    if len(expected) != len(records):
        raise ValueError("Duplicate alleles in VCF")
    results = {}
    for obj in vep_objects:
        fields = obj.get("input", "").split()
        if len(fields) < 5:
            raise ValueError("VEP JSON missing input VCF fields")
        key = (fields[0].removeprefix("chr"), int(fields[1]), fields[3].upper(), fields[4].upper())
        if key not in expected or key in results:
            raise ValueError(f"Unexpected or duplicate VEP annotation: {key}")
        if obj.get("assembly_name") not in (None, "GRCh38"):
            raise ValueError("VEP annotation assembly is not GRCh38")
        results[key] = obj
    if set(results) != set(expected):
        raise ValueError(f"VEP returned {len(results)}/{len(expected)} requested alleles")
    return [(key, summarize(obj, key[3]), obj) for key, obj in results.items()]


def summarize(obj, alt):
    transcripts = obj.get("transcript_consequences") or []
    # Keep all transcript annotations in raw JSON; one display transcript only.
    transcripts = sorted(transcripts, key=lambda t: (not bool(t.get("mane_select")),
                                                       not bool(t.get("canonical")),
                                                       t.get("transcript_id", "")))
    top = transcripts[0] if transcripts else {}
    frequencies = []
    for hit in obj.get("colocated_variants") or []:
        # A frequency is only meaningful for the matching alternate allele.
        per_allele = (hit.get("frequencies") or {}).get(alt) or {}
        value = per_allele.get("gnomadg")
        if value is not None:
            value = float(value)
            if not 0 <= value <= 1:
                raise ValueError("Invalid gnomAD genome AF")
            frequencies.append(value)
    if len(set(frequencies)) > 1:
        raise ValueError("Conflicting gnomAD frequencies for one allele")
    return {
        "gene_id": top.get("gene_id"),
        "gene_symbol": top.get("gene_symbol"),
        "transcript_id": top.get("transcript_id"),
        "consequence": obj.get("most_severe_consequence"),
        "gnomadg_af": frequencies[0] if frequencies else None,
    }


def load_jsonl(path):
    with open(path, encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip() and not line.startswith("#")]


def parse_hpo_gene_associations(path, gene_filter=None, limit=20000):
    """Read official genes_to_phenotype.txt with commented or plain header."""
    found = set()
    required = {"gene_symbol", "hpo_id", "hpo_name", "disease_id"}
    with open(path, encoding="utf-8-sig") as handle:
        header = None
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                maybe = line.lstrip("#").strip().split("\t")
                if required.issubset(maybe):
                    header = maybe
                continue
            if not header:
                maybe = line.rstrip("\r\n").split("\t")
                if required.issubset(maybe):
                    header = maybe
                    continue
                raise ValueError("HPO file missing expected tab-separated header")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                raise ValueError("Malformed HPO association")
            row = dict(zip(header, fields))
            gene = row["gene_symbol"]
            if gene_filter is not None and gene not in gene_filter:
                continue
            if not row["hpo_id"].startswith("HP:") or not row["disease_id"]:
                raise ValueError("Invalid HPO or disease identifier")
            found.add((gene, row["hpo_id"], row["hpo_name"], row["disease_id"]))
            if len(found) > limit:
                raise ValueError("HPO import limit exceeded")
    return sorted(found)


def load_clinvar_example(path):
    """Read a small, sourced condition assertion snapshot (not a live ClinVar feed)."""
    items = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(items, list) or not items or len(items) > 20:
        raise ValueError("Expected 1–20 ClinVar example assertions")
    required = {"assembly", "chrom", "pos", "ref", "alt", "gene_symbol", "condition",
                "disease_id", "classification", "review_status", "record_accession",
                "source_url", "snapshot_date"}
    for row in items:
        if not isinstance(row, dict) or set(row) != required or row["assembly"] != "GRCh38":
            raise ValueError("Invalid ClinVar assertion fields or assembly")
        if not isinstance(row["pos"], int) or row["pos"] < 1 or not all(
                isinstance(row[k], str) and row[k] for k in required - {"pos"}):
            raise ValueError("Invalid ClinVar assertion value")
        if not re.fullmatch(r"[ACGT]+", row["ref"]) or not re.fullmatch(r"[ACGT]+", row["alt"]):
            raise ValueError("Invalid ClinVar alleles")
        if not re.fullmatch(r"RCV\d+\.\d+", row["record_accession"]):
            raise ValueError("Invalid condition record accession")
        if not re.fullmatch(r"https://www\.ncbi\.nlm\.nih\.gov/clinvar/RCV\d+/", row["source_url"]):
            raise ValueError("ClinVar source must be an NCBI condition record")
        date.fromisoformat(row["snapshot_date"])
    return items
