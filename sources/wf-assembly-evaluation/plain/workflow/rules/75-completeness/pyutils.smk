import pathlib
import pandas


def create_mockup(file_path, caller):

    mock_name = pathlib.Path(file_path).name
    mock_path = DIR_LOCAL_REF.joinpath("mocks", mock_name)
    mock_path.parent.mkdir(exist_ok=True, parents=True)
    with open(mock_path, "w") as mock:
        _ = mock.write(f"{caller}\n{get_timestamp()}\n")
    return mock_path


def load_karyotype_estimate(karyo_est_file, sample, asm_unit):

    df = pandas.read_csv(karyo_est_file, sep="\t", header=0)
    select_unit = df["asm_unit"] == asm_unit
    assert select_unit.sum() == 1
    karyotype = df.loc[select_unit, "karyotype"].values[0]
    assert karyotype in ["male", "female", "any"]
    return karyotype


def get_karyotype_est_file(sample, refgenome):

    # look up estimated karyotype/sex for this assembly unit
    ### IMPORTANT
    # This path - see module 50-postprocess::estimate_asm_unit_karyotype
    ###
    fmt = {"sample": sample, "refgenome": refgenome}
    karyo_est_file = rules.estimate_asm_unit_karyotype.output.karyo_est
    karyo_est_file = str(karyo_est_file).format(**fmt)
    return karyo_est_file


def get_gene_model(gene_model, sample, refgenome, sex_or_unit):

    this_func = "75-completeness::pyutils::get_gene_model"
    temp_use_mockup = False
    if sex_or_unit in ["male", "female", "any"]:
        matched_model = config["gene_models"][gene_model][sex_or_unit]
    elif sex_or_unit not in ASSEMBLY_UNITS_SEX_SPECIFIC:
        # other parts of the assembly that do not include
        # any of the sex chromosomes
        matched_model = config["gene_models"][gene_model]["any"]
    else:
        assert sample in SAMPLES
        assert sex_or_unit in ASSEMBLY_UNITS_SEX_SPECIFIC
        karyo_est_file = get_karyotype_est_file(sample, refgenome)
        if not pathlib.Path(karyo_est_file).is_file():
            # file has not been created yet - return mock
            # to make this an explicit failure
            matched_model = f"{sample}-{sex_or_unit}-{refgenome}-no-karyotype.mock"
            temp_use_mockup = True
        else:
            karyotype_est = load_karyotype_estimate(karyo_est_file, sample, sex_or_unit)
            matched_model = config["gene_models"][gene_model][karyotype_est]
    model_path = DIR_GLOBAL_REF.joinpath(matched_model)
    if temp_use_mockup:
        model_path = create_mockup(model_path, this_func)
    return model_path


def get_reference_gene_model_alignment(gene_model, sample, refgenome, asm_unit):

    this_func = "75-completeness::pyutils::get_reference_gene_model_alignment"
    temp_use_mockup = False
    karyo_est_file = get_karyotype_est_file(sample, refgenome)
    if not pathlib.Path(karyo_est_file).is_file() and asm_unit in ASSEMBLY_UNITS_SEX_SPECIFIC:
        # file has not been created yet - return mock
        # to make this an explicit failure
        matched_alignment = f"{sample}-{asm_unit}-{refgenome}-no-karyotype.mock"
        temp_use_mockup = True
    elif asm_unit not in ASSEMBLY_UNITS_SEX_SPECIFIC:
        # some other assembly unit / not sex-specific - just use 'any' here
        match_alignment_fmt = {"refgenome": refgenome, "genemodel": gene_model}
        matched_alignment = str(rules.ref_completeness_genemodel_any.output.paf)
        matched_alignment = matched_alignment.format(**match_alignment_fmt)
    else:
        karyotype_est = load_karyotype_estimate(karyo_est_file, sample, asm_unit)
        match_alignment_fmt = {"refgenome": refgenome, "genemodel": gene_model}
        if karyotype_est == "female":
            matched_alignment = str(rules.ref_completeness_genemodel_female.output.paf)
        elif karyotype_est == "male":
            matched_alignment = str(rules.ref_completeness_genemodel_male.output.paf)
        else:
            matched_alignment = str(rules.ref_completeness_genemodel_any.output.paf)
        matched_alignment = matched_alignment.format(**match_alignment_fmt)
    if temp_use_mockup:
        matched_alignment = create_mockup(matched_alignment, this_func)
    return matched_alignment


def get_reference_genome(refgenome, sex="any"):

    try:
        matched_ref = config["refgenomes"][refgenome][sex]
    except KeyError:
        matched_ref = config["refgenomes"][refgenome]["any"]
    ref_path = DIR_GLOBAL_REF.joinpath(matched_ref)
    return ref_path


def get_contig_to_reference_norm_paf(wildcards):

    _this_func = "75-completeness::pyutils::get_contig_to_reference_norm_paf"

    assert hasattr(wildcards, "aln_type"), f"Invalid wildcards - no 'aln_type': {_this_func} / {wildcards}"

    if wildcards.aln_type == "blevel":
        norm_paf = rules.normalize_minimap_assembly_to_reference_align_paf.output.tsv
    elif wildcards.aln_type == "coarse":
        norm_paf = rules.normalize_mashmap_assembly_to_reference_align_paf.output.tsv
    else:
        err_msg = (
            f"ERROR in function: {_this_func}\n"
            "Expected value 'blevel' or 'coarse' for wildcard 'aln_type'\n"
            f"Received: {wildcards.aln_type} / {wildcards}"
        )
        logerr(err_msg)
        raise ValueError(err_msg)

    return norm_paf
