
def get_reference_assembly(sample, reference):

    sample_sex = SAMPLE_INFOS[sample]["sex"]
    try:
        matched_ref = config["refgenomes"][reference][sample_sex]
    except KeyError:
        matched_ref = config["refgenomes"][reference]["any"]
    ref_path = DIR_GLOBAL_REF.joinpath(matched_ref)
    return ref_path
