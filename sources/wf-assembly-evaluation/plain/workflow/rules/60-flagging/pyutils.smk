

def flatten_window_readcov_columns(column_level_names, columns):

    aln_types = {
        0: "SPL",
        1: "PRI",
        2: "SEC"
    }
    assert column_level_names[0] == "read_type"
    assert column_level_names[1] == "aln_type"
    assert column_level_names[2] == "mapq"
    assert column_level_names[3] == "value"

    flattened_columns = []
    for column_tuple in columns:
        if column_tuple[3] in ["name", "start", "end", "contig"]:
            flattened_columns.append(column_tuple[3])
        else:
            flat_column = (
                f"{column_tuple[0]}_"
                f"ALN-{aln_types[column_tuple[1]]}_"
                f"MAPQ-{column_tuple[2]}_"
                f"{column_tuple[3]}"
            )
            flattened_columns.append(flat_column)

    return flattened_columns


def load_reference_centromere_annotation(ref_name):

    try:
        annotation = config["refgenomes"][ref_name]["censat"]
    except KeyError:
        logerr(f"No centromere annotation found in config for: {ref_name}")
        raise ValueError(f"No centromere annotation found in config: {ref_name}")
    file_path = DIR_GLOBAL_REF.joinpath(annotation)
    return file_path
