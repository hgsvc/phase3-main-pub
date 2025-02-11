import pandas


def load_assembly_unit_karyotypes(file_path):

    df = pandas.read_csv(file_path, sep="\t", comment="#", header=0)

    lut = dict()
    for sample, au_infos in df.groupby("sample"):
        sex_hap1 = au_infos["karyotype"].iloc[0]
        sex_hap2 = au_infos["karyotype"].iloc[1]
        assert au_infos["asm_unit"].iloc[0] == "asm-hap1"
        # this is a manual fix for sample HG00732 that
        # lacks a decent chrX assembly due to cell line
        # issues
        if sex_hap1 == sex_hap2:
            sex_label = sex_hap1
        elif sample.startswith("HG00732"):
            sex_label = "female"
        else:
            sex_label = f"{sex_hap1},{sex_hap2}"
        lut[sample] = sex_label
    return lut
