

def compute_karyotype_score(chrom_assign_table, karyo_chroms):
    """
    Take only alignments into account that cover at least 10
    percent of the query sequence length. Then, sum up over
    matching bases over all alignments
    """
    select_chroms = chrom_assign_table["target_name"].isin(karyo_chroms)
    select_cov = chrom_assign_table["match_cov_query"] > 0.1
    selector = select_chroms & select_cov
    subset = chrom_assign_table.loc[selector, :].copy()
    if subset.empty:
        karyo_score = 0
        karyo_fragmentation = 0
        karyo_contigs = 0
    else:
        karyo_score = int(subset["align_matching"].sum())
        karyo_fragmentation = int(subset["align_num"].sum())
        karyo_contigs = int(subset["query_name"].nunique())
    return karyo_score, karyo_fragmentation, karyo_contigs

