import unittest
from mhca import annotate_haplotype 
from pathlib import Path
import mhca
from argparse import Namespace
import filecmp
import gzip

class TestAnnotate(unittest.TestCase):
    def setUp(self):
        p = Path(mhca.__file__)
        testfolder = p.parent.parent / 'tests' 
        zippath = testfolder / 'test1.fasta.gz'
        self.fastapath = testfolder / 'test1.fasta'
        with gzip.open(zippath, 'rb') as fasta_zipped, open(self.fastapath, 'wb') as fasta_plain:
            file_content = fasta_zipped.read()
            fasta_plain.write(file_content)
        self.outfolder = testfolder / 'tmp_out'
        self.outfolder.mkdir()
        self.outfile_template = testfolder / 'test1_template.gff'
        self.outfile_test = self.outfolder / 'test1.gff'

    def tearDown(self):
        self.fastapath.unlink()
        for result_file in self.outfolder.iterdir():
            result_file.unlink()
        self.outfolder.rmdir()
    
    def test_parse_cigar(self):
        result = annotate_haplotype.parse_cigar("3M2D10M", "+", "AA|CCGGTT|AACCGG|TTAACCGG|TT", 3)
        self.assertEqual(result, [(0, 7), (8, 13)])

    def test_without_rs(self):
        args = Namespace(   haplotype = str(self.fastapath), \
                            manual_corrections = None, \
                            locus_tag_prefix = None, \
                            imgt_folder = None, \
                            output_folder = self.outfolder, \
                            skip_mapping = False, \
                            skip_imgt = False, \
                            skip_rsg = False, \
                            skip_rs = True, \
                            refseq_full_fasta = None, \
                            refseqgene_full_fasta = None)
        annotate_haplotype.main(args)
        self.assertTrue(self.outfile_test.is_file())
        self.assertTrue(filecmp.cmp(self.outfile_test,self.outfile_template,shallow=True))
        

if __name__ == '__main__':
    unittest.main()
