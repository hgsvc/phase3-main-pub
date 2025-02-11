#!/usr/bin/env python

import yaml

addSample = True
SM_dict = {}


while addSample == True:
    outfile = input("What do you want to call the config file: ")
    SM = input("Input sample name: ")
    SM_dict[SM] = {}

    HiFi_FOFN = input("Input HiFi_FOFN if none press enter: ")
    ONT_FOFN = input("Input ONT_FOFN if none press enter: ")
    SM_dict[SM]["fofns"] = {}
    SM_dict[SM]["type_map"] = {}

    if HiFi_FOFN != "":
        type_map_hifi = input(
            "map HiFi reads with minimap2 or winnowmap? (type minimap2 or winnowmap): "
        )
        SM_dict[SM]["fofns"]["HiFi"] = HiFi_FOFN
        SM_dict[SM]["type_map"]["HiFi"] = type_map_hifi
    if ONT_FOFN != "":
        type_map_ont = input(
            "map ONT reads with minimap2 or winnowmap? (type minimap2 or winnowmap): "
        )
        SM_dict[SM]["fofns"]["ONT"] = ONT_FOFN
        SM_dict[SM]["type_map"]["ONT"] = type_map_ont

    SM_dict[SM]["regions"] = {}
    addRegion = True
    while addRegion == True:
        region = input("name region, if none press enter: ")
        if region == '':
            break
        bed = input("Path to bed file of contigs in named region: ")
        SM_dict[SM]["regions"][region] = bed
        nextRegion = input("Add another region? (y)es, (n)o: ")
        if nextRegion == "n":
            addRegion = False
    repeatMask = input("Run Repeat Mask? (y)es, (n)o: ")
    if repeatMask == "y":
        SM_dict[SM]["repeat_mask"] = True
        species = input("Species for Repeat Mask?: ")
        SM_dict[SM]["species"] = species
    else:
        SM_dict[SM]["repeat_mask"] = False
        SM_dict[SM]["species"] = "human"

    genBed = input("Generate bed of potential collapses/misassemblies (y)es, (n)o: ")
    if genBed == "y":
        SM_dict[SM]["genBed"] = True
    else:
        SM_dict[SM]["genBed"] = False

    nucopt = input("flags for nucFreq? i.e. -y 100 if none press enter: ")
    if nucopt == "":
        SM_dict[SM]["nuc_opts"] = "-y 100"
    else:
        SM_dict[SM]["nuc_opts"] = nucopt
    hap1 = input("Path to hap1 of assembly: ")
    hap2 = input("Path to hap2 of assembly, if already combined, press enter: ")
    if hap2 == "":
        hap2 = "empty.fa"
        open(hap2, "a+").close()
    SM_dict[SM]["asm_h1"] = hap1
    SM_dict[SM]["asm_h2"] = hap2
    extra_tigs = input("If you have additional contigs which should be considered in the alignment, you can provide them with an fofn, enter path to fofn or press enter for none: ")
    if extra_tigs != "":
        SM_dict[SM]["extra_tigs"] = extra_tigs
    nextSample = input("add another sample? (y)es, (n)o: ")
    if nextSample == "n":
        addSample = False

with open(outfile, "w+") as build_yaml:
    yaml.dump(SM_dict, build_yaml, default_flow_style=False)
