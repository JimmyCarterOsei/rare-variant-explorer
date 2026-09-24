import tempfile
import unittest
from pathlib import Path

from explorer.processing import (annotations_for_records, load_clinvar_example, parse_hpo_gene_associations,
                                 parse_vcf, summarize)


class ProcessingTests(unittest.TestCase):
    def write(self, content):
        temp = tempfile.NamedTemporaryFile(mode="w", delete=False)
        temp.write(content)
        temp.close()
        self.addCleanup(Path(temp.name).unlink)
        return temp.name

    def test_vcf_requires_explicit_assembly_and_splits_alleles(self):
        vcf = "##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t20\t.\tA\tC,G\t.\t.\t.\n"
        records = parse_vcf(self.write(vcf))
        self.assertEqual([r[3] for r in records], ["C", "G"])
        self.assertEqual(records[1][4].split()[4], "G")
        with self.assertRaisesRegex(ValueError, "GRCh38"):
            parse_vcf(self.write(vcf.replace("##reference=GRCh38", "##reference=GRCh37")))

    def test_vep_requires_exact_allele_match(self):
        rows = [("1", 20, "A", "C", "1\t20\t.\tA\tC\t.\t.\t.")]
        response = {"input": rows[0][4], "most_severe_consequence": "missense_variant"}
        self.assertEqual(len(annotations_for_records(rows, [response])), 1)
        with self.assertRaisesRegex(ValueError, "Unexpected"):
            annotations_for_records(rows, [{**response, "input": response["input"].replace("\tC\t", "\tG\t")}])
        with self.assertRaisesRegex(ValueError, "1/2"):
            annotations_for_records(rows + [("1", 21, "A", "C", "1 21 . A C . . .")], [response])

    def test_frequency_is_allele_specific_and_missing_stays_null(self):
        obj = {"transcript_consequences": [{"gene_symbol": "ABC", "canonical": 1},
                                             {"gene_symbol": "XYZ", "mane_select": "NM_1"}],
               "colocated_variants": [{"frequencies": {"G": {"gnomadg": .31},
                                                        "C": {"gnomadg": 0.0}}}]}
        result = summarize(obj, "C")
        self.assertEqual(result["gene_symbol"], "XYZ")
        self.assertEqual(result["gnomadg_af"], 0.0)
        self.assertIsNone(summarize(obj, "T")["gnomadg_af"])
        with self.assertRaisesRegex(ValueError, "Conflicting"):
            summarize({"colocated_variants": [{"frequencies": {"C": {"gnomadg": .1}}},
                                             {"frequencies": {"C": {"gnomadg": .2}}}]}, "C")

    def test_hpo_import_keeps_disease_provenance(self):
        path = self.write("#ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n"
                          "16\tAARS1\tHP:0002460\tDistal muscle weakness\t15/15\tOMIM:613287\n"
                          "16\tAARS1\tHP:0002460\tDistal muscle weakness\t1/2\tORPHA:123\n")
        rows = parse_hpo_gene_associations(path, {"AARS1"})
        self.assertEqual(len(rows), 2)
        self.assertEqual({r[3] for r in rows}, {"OMIM:613287", "ORPHA:123"})
        self.assertEqual(parse_hpo_gene_associations(path, {"OTHER"}), [])

    def test_hpo_plain_header_after_comments(self):
        path = self.write("#downloaded from HPO\n"
                          "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n"
                          "1080\tCFTR\tHP:0000007\tAutosomal recessive inheritance\t-\tOMIM:219700\n")
        self.assertEqual(parse_hpo_gene_associations(path, {"CFTR"}),
                         [("CFTR", "HP:0000007", "Autosomal recessive inheritance", "OMIM:219700")])

    def test_cftr_vcf_and_condition_snapshot_match(self):
        root = Path(__file__).resolve().parent.parent
        record = parse_vcf(root / "data/cftr_f508del.vcf")[0]
        assertion = load_clinvar_example(root / "data/clinvar_example.json")[0]
        self.assertEqual(record[:4], (assertion["chrom"], assertion["pos"],
                                      assertion["ref"], assertion["alt"]))
        self.assertEqual(assertion["disease_id"], "OMIM:219700")
        self.assertTrue(assertion["record_accession"].startswith("RCV"))

    def test_cftr_g551d_vcf_and_clinvar_snapshot_match(self):
        root = Path(__file__).resolve().parent.parent
        record = parse_vcf(root / "data/cftr_g551d.vcf")[0]
        assertion = load_clinvar_example(root / "data/clinvar_g551d.json")[0]
        self.assertEqual(record[:4], (assertion["chrom"], assertion["pos"],
                                      assertion["ref"], assertion["alt"]))
        self.assertEqual(assertion["record_accession"], "RCV000007540.61")


if __name__ == "__main__":
    unittest.main()
