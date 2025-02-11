
ACROCENTRIC_CHROMOSOMES = ["chr13", "chr14", "chr15", "chr21", "chr22"]
ACROS = ACROCENTRIC_CHROMOSOMES


# chromosome extraction - specific settings
CHROM_Y_SELECT_MOTIFS = [
    "DYZ18_Yq",
    "DYZ1_Yq",
    "DYZ2_Yq",
    "DYZ3-sec_Ycentro"
]

CHROM_Y_ALL_MOTIFS = [
    "DYZ18_Yq",
    "DYZ1_Yq",
    "DYZ2_Yq",
    "DYZ3-sec_Ycentro",
    "DYZ19_Yq",
    "DYZ3-prim_Ycentro",
    "TSPY",
    "Yqhet_2k7bp",
    "Yqhet_3k1bp",
]

HMMER_MOTIF_PARAMS = {
    "DYZ18_Yq": {
        "score_t": 2100,
        "evalue_t": '1.60E-150'
    },
    "DYZ1_Yq": {
        "score_t": 2500,
        "evalue_t": '1.60E-150'
    },
    "DYZ2_Yq": {
        "score_t": 1700,
        "evalue_t": '1.60E-150',
        "scale_cpu": 2,
        "scale_time": 4
    },
    "DYZ3-sec_Ycentro": {
        "score_t": 1700,
        "evalue_t": '1.60E-150'
    },
    "TSPY": {
        "score_t": 1000,
        "evalue_t": '1.60E-200',
        "scale_cpu": 4,
        "scale_time": 47,
        "scale_mem": 7
    },
    "DYZ19_Yq": {
        "evalue_t": '1.60E-15'
    },
    "DYZ3-prim_Ycentro": {
        "score_t": 90,
        "evalue_t": '1.60E-15'
    },
    "Yqhet_2k7bp": {
        "score_t": 1500,
        "evalue_t": '1.60E-200'
    },
    "Yqhet_3k1bp": {
        "score_t": 1500,
        "evalue_t": '1.60E-200'
    }
}

HMMER_MOTIF_NAMES = sorted(HMMER_MOTIF_PARAMS.keys())
