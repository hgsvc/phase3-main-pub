import pandas

SAMPLE_TABLE = pandas.read_csv(
    config["samples"],
    sep="\t",
    comment="#",
    header=0
)

SAMPLES = sorted(SAMPLE_TABLE["sample"].unique())

PLAIN_SAMPLES = [s.split(".")[0] for s in SAMPLES]
assert len(PLAIN_SAMPLES) == len(SAMPLES)

SAMPLE_SEX = None
MALE_SAMPLES = None
FEMALE_SAMPLES = None
if "sex" in SAMPLE_TABLE:
    SAMPLE_SEX = dict()
    for row in SAMPLE_TABLE.itertuples():
        SAMPLE_SEX[row.sample] = row.sex
    MALE_SAMPLES = [s for s in SAMPLES if SAMPLE_SEX[s] in ["male", "m"]]
    FEMALE_SAMPLES = [s for s in SAMPLES if SAMPLE_SEX[s] in ["female", "f"]]

if all("vrk-ps" in s for s in SAMPLES):

    ASSEMBLER = "verkko"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2", "unassigned"]

elif all(("CEPH" in s and "NA" in s) for s in SAMPLES):

    # special setting for Verkko pedigree samples
    ASSEMBLER = "verkko"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2", "unassigned"]

elif all("CEPH" in s for s in SAMPLES):

    # special setting for hifiasm pedigree samples
    ASSEMBLER = "hifiasm"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2"]

elif all("hsm-ps" in s for s in SAMPLES):

    ASSEMBLER = "hifiasm"
    MAIN_ASSEMBLY_UNITS = ["hap1", "hap2"]

else:
    raise ValueError("Can process either Verkko or hifiasm assemblies, not both.")


def get_asm_unit_fasta_files(wildcards):
    """Used in:
    module extract_chrom::10_preselect_tigs.smk (per sample)
    module regions::hla::extract.smk (all samples)
    """
    assert ASSEMBLER in ["verkko", "hifiasm"]
    if ASSEMBLER == "verkko":
        if hasattr(wildcards, "sample"):
            # used in chrY / extract
            return expand(
                WORKDIR_EVAL.joinpath(
                    "results/assemblies/", "{{sample}}",
                    "{{sample}}.asm-{asm_unit}.fasta.gz"
                ),
                asm_unit=MAIN_ASSEMBLY_UNITS
            )
        else:
            # used in HLA / extract
            return expand(
                WORKDIR_EVAL.joinpath(
                    "results/assemblies", "{sample}",
                    "{sample}.asm-{asm_unit}.fasta.gz"
                ),
                sample=SAMPLES,
                asm_unit=MAIN_ASSEMBLY_UNITS
            )
    else:
        fasta_columns = [f"asm_{au}" for au in MAIN_ASSEMBLY_UNITS]
        if hasattr(wildcards, "sample"):
            fasta_files = SAMPLE_TABLE.loc[SAMPLE_TABLE["sample"] == wildcards.sample, fasta_columns].values[0]
            fasta_files = [pathlib.Path(f).resolve(strict=True) for f in fasta_files]
        else:
            sub = SAMPLE_TABLE[fasta_columns].values
            fasta_files = []
            for pair in sub:
                f1, f2 = pair
                fasta_files.append(
                    pathlib.Path(f1).resolve(strict=True)
                )
                fasta_files.append(
                    pathlib.Path(f2).resolve(strict=True)
                )
        return fasta_files
