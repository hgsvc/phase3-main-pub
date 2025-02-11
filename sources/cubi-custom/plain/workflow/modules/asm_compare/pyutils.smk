

def to_int(number):

    try:
        integer = int(number)
    except ValueError:
        suffix = number[-1]
        numeric = int(number[:-1])
        if suffix.lower() == "k":
            factor = 1e3
        elif suffix.lower() == "m":
            factor = 1e6
        else:
            raise ValueError(number)
        integer = int(numeric * factor)
    return integer


def parse_precision_param(settings, get):

    try:
        setting, threshold = settings.split("-")
        threshold = to_int(threshold)
    except ValueError:
        assert settings == "exact"
        setting = "exact"
        threshold = 50

    if get == "setting":
        return setting
    elif get == "threshold":
        return threshold
    else:
        raise ValueError(get)


def select_genome_size_file(wildcards):

    if hasattr(wildcards, "qry_to_trg"):
        if wildcards.qry_to_trg == "vrk-to-hsm":
            fai = rules.align_verkko_to_hifiasm.input.hsm_idx
        elif wildcards.qry_to_trg == "hsm-to-vrk":
            fai = rules.align_hifiasm_to_verkko.input.vrk_idx
        else:
            raise ValueError(wildcards)
    elif hasattr(wildcards, "child"):
        fai = rules.align_verkko_parent_to_child.input.chl_idx
    else:
        raise ValueError(wildcards)
    return fai
