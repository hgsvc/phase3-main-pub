from mhca import annotate_haplotype, check_CDS, refseq2fullfasta, update_feature_table
from sys import argv, exit
from argparse import ArgumentParser

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help', dest="subparser_name")

    anno_parser = subparsers.add_parser("annotate")
    anno_parser.add_argument("haplotype", help="Input Haplotype in fasta format.")
    anno_parser.add_argument("--skip_imgt", default=False, action="store_true", help="Use this if you don't want to use IMGT as a data source.")
    anno_parser.add_argument("--skip_rsg", default=False, action="store_true", help="Use this if you don't want to use RefSeqGene as a data source.")
    anno_parser.add_argument("--skip_rs", default=False, action="store_true", help="Use this if you don't want to use RefSeq as a data source.")
    
    anno_parser.add_argument("--skip_mapping", default=False, action="store_true", help="If you already have used minimap2 and want to skip this step.")
    anno_parser.add_argument("--locus_tag_prefix", help="Comma separated file containing the locus_tag_prefix for the haplotype")
    anno_parser.add_argument("--manual_corrections", help="Comma separated file with manual corrections")
    anno_parser.add_argument("--imgt_folder", help="Folder which holds the start to stop codon fastas and the full fastas with exon boundary information.")
    anno_parser.add_argument("--refseqgene_full_fasta", help="Fasta file with exon boundary information (boundaries denoted as: '|'). This can be generated from the RefSeqGene genbank file and the genbank2fullfasta method of this package. If this parameter is omitted, a precompiled version of this is used.")
    anno_parser.add_argument("--refseq_full_fasta", help="Fasta file with exon boundary information (boundaries denoted as: '|'). This can be generated from a RefSeq dataset. For details have a look at the refseq2fullfasta script of this package. If this parameter is omitted, a precompiled version of this is used.")
    anno_parser.add_argument("output_folder")

    CDS_parser = subparsers.add_parser("check_CDS")
    CDS_parser.add_argument("haplotype", help="Input Haplotype in fasta format.")
    CDS_parser.add_argument("annotation_gff")

    refseq_parser = subparsers.add_parser("refseq2fullfasta", help="Needs samtools faidx")
    refseq_parser.add_argument("reference", help="Genome reference file. Usually hg38 in fasta format")
    refseq_parser.add_argument("refseq_genes", help="Refseq tab-separated file")
    refseq_parser.add_argument("outfile", help="Output fasta file.")

    update_feature_parser = subparsers.add_parser("update_feature_table")
    update_feature_parser.add_argument("old_feature_table", help="Old feature table file.")
    update_feature_parser.add_argument("record_id", help="ID of record that will be updated")
    update_feature_parser.add_argument("update_feature_table", help="Feature table file with new features.")
    update_feature_parser.add_argument("new_feature_table", help="New updated feature table file.")

    args = parser.parse_args()
    if args.subparser_name == "annotate":
        return annotate_haplotype.main( args )
    elif args.subparser_name == "check_CDS":
        return check_CDS.main( args )
    elif args.subparser_name == "refseq2fullfasta":
        return refseq2fullfasta.main( args )
    elif args.subparser_name == "update_feature_table":
        return update_feature_table.main( args )
    else:
        parser.print_help()
