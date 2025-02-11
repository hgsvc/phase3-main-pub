#!/usr/bin/env python3

import argparse as argp
import collections as col
import functools as fnt
import hashlib as hl
import pathlib as pl

import pandas as pd
import networkx as nx


__prog__ = "merge_verkko_infos.py"
__developer__ = "Peter Ebert"
__license__ = "MIT"
__version__ = "(prototype)"


# TODO: make dataclass
CONTIG_UNITIG_PIECES = {
    "unitig": None,
    "contig": None,
    "hap": "H0",
    "is_used": "undefined",
    "layout_assigned": "undefined",
    "contig_pieces": None,
    "pieces_path": None,
    "total_gap_length": 0,
    "default_gaps": 0,
    "min_gaps": 0,
    "ambiguous": 0
}


def parse_command_line():

    parser = argp.ArgumentParser(prog=__prog__)
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=f"{__prog__} v{__version__}"
    )

    parser.add_argument(
        "--graph",
        "-g",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="graph",
        help="Path to Verkko output graph (GFA; assembly.homopolymer-compressed.gfa)."
    )
    parser.add_argument(
        "--fasta",
        "-f",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None,
        dest="fasta",
        help="Path to assembly FASTA to add contig length in uncompressed bp (assembly.fasta)."
    )
    parser.add_argument(
        "--scf-layout",
        "-sl",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="layout",
        help="Path to Verkko layout file (6-layoutContigs/unitig-popped.layout.scfmap)."
    )
    parser.add_argument(
        "--hifi-node-cov",
        "-hnc",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="hifi_cov",
        help="Path to Verkko node HIFI coverage file (assembly.hifi-coverage.csv)."
    )
    parser.add_argument(
        "--ont-node-cov",
        "-onc",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="ont_cov",
        help="Path to Verkko node ONT coverage file (assembly.ont-coverage.csv)."
    )
    parser.add_argument(
        "--rukki-colors",
        "-colors",
        "-rc",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None,
        dest="rukki_colors",
        help="Optional: path to Rukki node assignment/color output file (6-rukki/unitig-popped-unitig-normal-connected-tip.colors.csv)"
    )
    parser.add_argument(
        "--rukki-final",
        "-final",
        "-rf",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None,
        dest="rukki_final",
        help="Optional: path to intermediate Rukki node assignment/color output file (6-rukki/out_final_ann.csv)"
    )
    parser.add_argument(
        "--rukki-paths",
        "-paths",
        "-rp",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None,
        dest="rukki_paths",
        help="Optional: path to Rukki paths output file (6-rukki/unitig-popped-unitig-normal-connected-tip.paths.tsv)"

    )
    parser.add_argument(
        "--out-table",
        "-table",
        "-ot",
        type=lambda x: pl.Path(x).resolve(strict=False),
        required=True,
        dest="table",
        help="Path to output table listing all unitigs with merged infos (TSV)."
    )
    parser.add_argument(
        "--out-path-ids",
        "-ids",
        "-opi",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None,
        dest="path_ids",
        help="Path to output table to dump rukki paths plus path IDs (TSV)."
    )

    args = parser.parse_args()

    return args


def read_coverage_file(file_path, data_source):

    cov_table = pd.read_csv(
        file_path,
        header=0,
        sep="\t"
    )
    cov_table[f"{data_source}_cov"] = cov_table["coverage"].round(2)
    cov_table.drop("coverage", axis=1, inplace=True)
    return cov_table


def parse_layout_unitig_info(unitig_info):
    """Helper function to deal with different
    information in layout file depending on
    phasing state of the assembly.

    Args:
        unitig_info (_type_): _description_

    Returns:
        _type_: _description_
    """
    infos = dict(CONTIG_UNITIG_PIECES)
    if unitig_info.startswith("utig"):
        # not phased/scaffolded
        infos["unitig"] = unitig_info
    else:
        parts = unitig_info.split("_")
        assert len(parts) == 3, f"Cannot process layout: {parts}"
        assert parts[2].startswith("utig"), f"Cannot process layout: {parts}"
        infos["unitig"] = parts[2]
        if parts[1] == "from":
            infos["is_used"] = "yes"
        if parts[0] == "mat":
            infos["layout_assigned"] = "maternal"
            infos["hap"] = "H1"
        elif parts[0] == "pat":
            infos["layout_assigned"] = "paternal"
            infos["hap"] = "H2"
        else:
            assert parts[0] == "na"
    return infos


def read_layout_file(file_path):

    current_utg = None
    tig_mapping = dict()
    with open(file_path, "r") as listing:
        for ln, line in enumerate(listing, start=1):
            if not line.strip():
                continue
            if line.startswith("path"):
                _, contig, tig_desc = line.strip().split()
                tig_infos = parse_layout_unitig_info(tig_desc)
                tig_infos["contig"] = contig
                current_utg = tig_infos["unitig"]
            elif line.startswith("piece"):
                piece = line.strip()
                if current_utg in tig_mapping:
                    first_entry = tig_mapping[current_utg]
                    first_entry["ambiguous"] = 1
                    first_entry["contig"] += f"|{tig_infos['contig']}"
                    first_entry["hap"] = "H0"
                    first_entry["layout_assigned"] = "undefined"
                    first_entry["contig_pieces"] = first_entry["contig_pieces"].split(",")
                    first_entry["pieces_path"] = first_entry["pieces_path"].split("->-")
                    tig_infos = first_entry
                #assert current_utg not in tig_mapping, f"Duplicate: {current_utg}"
                assert current_utg is not None
                try:
                    tig_infos["contig_pieces"].append(piece)
                    tig_infos["pieces_path"].append(piece)
                except AttributeError as ae:
                    tig_infos["contig_pieces"] = [piece]
                    tig_infos["pieces_path"] = [piece]
            elif line.startswith("[N"):
                gap_size = int(line.strip().strip("[]N"))
                tig_infos["pieces_path"].append(f"gap{gap_size}")
                if gap_size == 1000:
                    tig_infos["min_gaps"] += 1
                elif gap_size == 5000:
                    tig_infos["default_gaps"] += 1
                else:
                    pass
            elif line.startswith("end"):
                if tig_infos["ambiguous"] == 1:
                    tig_infos["contig_pieces"] = "|".join(tig_infos["contig_pieces"])
                    tig_infos["pieces_path"] = "|".join(tig_infos["pieces_path"])
                else:
                    tig_infos["contig_pieces"] = ",".join(tig_infos["contig_pieces"])
                    tig_infos["pieces_path"] = "->-".join(tig_infos["pieces_path"])
                tig_mapping[tig_infos["unitig"]] = dict(tig_infos)
                tig_infos = None
                current_utg = None
            else:
                raise ValueError(f"{ln}: {line.strip()}")
    # means file was ended with "end" line
    assert current_utg is None, "Incomplete parsing of layout file"
    return tig_mapping


def build_graph(graph_file):

    graph = nx.Graph()
    node_lengths = col.Counter()
    with open(graph_file, "r") as gfa:
        for ln, line in enumerate(gfa, start=1):
            if line.startswith("S"):
                columns = line.strip().split()
                node = columns[1]
                if columns[2] == "*":
                    ln_tag = [c for c in columns[3:] if c.startswith("LN")]
                    assert len(ln_tag) == 1
                    # entry of form: LN:i:1761044
                    node_length = int(ln_tag[0].split(":")[-1])
                else:
                    node_length = len(columns[2])
                node_lengths[node] = node_length
                graph.add_node(node)
            elif line.startswith("L"):
                columns = line.strip().split()
                node_a = columns[1]
                node_b = columns[3]
                graph.add_edge(node_a, node_b)
            else:
                raise ValueError(f"{ln}: {line.strip()}")
    return graph, node_lengths


def process_connected_components(graph, node_lengths, layout):

    unitig_records = []
    cc_records = []
    for cc_id, cc_nodes in enumerate(nx.connected_components(graph), start=1):
        cc = graph.subgraph(cc_nodes).copy()
        size_of_cc = len(cc.nodes)
        length_of_cc = sum(node_lengths[n] for n in cc.nodes)
        cc_record = {
            "cc_id": f"CC{cc_id}",
            "cc_size": size_of_cc,
            "cc_length_hpc": length_of_cc
        }
        cc_records.append(cc_record)
        for node in cc.nodes:
            # node = unitig
            unitig_record = dict(cc_record)
            unitig_record["node"] = node
            unitig_record["node_num"] = int(node.split("-")[-1])
            node_length = node_lengths[node]
            unitig_record["node_length_hpc"] = node_length
            try:
                layout_info = layout[node]
            except KeyError:
                assert node not in layout
                layout_info = dict(CONTIG_UNITIG_PIECES)
                layout_info["unitig"] = node
                layout_info["contig"] = "no-contig"
                layout_info["contig_pieces"] = "no-layout"
                layout_info["pieces_path"] = "no-layout"
            unitig_record.update(layout_info)
            unitig_records.append(unitig_record)

    utg_df = pd.DataFrame.from_records(unitig_records)
    assert not pd.isnull(utg_df).all(axis=0).all()

    return utg_df


def compute_rukki_path_ids(path_file):

    df = pd.read_csv(path_file, sep="\t", header=0)

    nodes_to_id = dict()
    path_ids = []
    for row in df.itertuples(index=False):
        path_id = hl.md5(row.path.encode("-utf-8")).hexdigest()
        path_ids.append(path_id)
        for utig in row.path.split(","):
            if "utig" not in utig:
                continue
            nodes_to_id[utig.strip("+-")] = path_id
    df["path_id"] = path_ids
    return df, nodes_to_id


def process_final_node_assignment_table(unitigs, assign_path):

    final_assign = pd.read_csv(assign_path, sep="\t", header=0)
    final_assign.drop(["length", "info"], axis=1, inplace=True)
    final_assign["node_class"] = final_assign["assignment"].str.lower()
    final_assign["class_color"] = final_assign["color"]
    final_assign.drop(["assignment", "color"], axis=1, inplace=True)
    unitigs = unitigs.merge(final_assign, on="node", how="outer")
    unitigs["node_class"].fillna("unassigned", inplace=True)
    unitigs["class_color"].fillna("#000000", inplace=True)

    return unitigs


def process_rukki_color_table(unitigs, color_path):

    rukki_colors = pd.read_csv(color_path, sep="\t", header=0)
    rukki_colors["rukki_color"] = rukki_colors["color"]
    rukki_colors["mat_marker"] = rukki_colors["mat"]
    rukki_colors["pat_marker"] = rukki_colors["pat"]
    rukki_colors["mat_div_pat"] = rukki_colors["mat:pat"]
    rukki_colors.drop(["color", "mat", "pat", "mat:pat"], axis=1, inplace=True)
    unitigs = unitigs.merge(rukki_colors, on="node", how="outer")
    unitigs["rukki_color"].fillna("#000000", inplace=True)
    unitigs["mat_marker"].fillna(0, inplace=True)
    unitigs["pat_marker"].fillna(0, inplace=True)
    unitigs["mat_div_pat"].fillna("0:0", inplace=True)

    return unitigs


def read_sequence_lengths(fasta_path):

    contig_lengths = col.Counter()
    current_contig = None
    contig_length = 0
    with open(fasta_path, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                if current_contig is not None:
                    contig_lengths[current_contig] = contig_length
                current_contig = line.strip()[1:]
                contig_length = 0
            elif line.strip():
                contig_length += len(line)
            else:
                continue
    contig_lengths[current_contig] = contig_length
    return contig_lengths


def assign_contig_length(contig_lengths, contig):
    """This function is necessary to deal with the
    ambiguous assignments that occur every now and then:

    haplotype1-0000484|haplotype2-0002728

    Args:
        contig (_type_): _description_
        contig_lengths (_type_): _description_
    """

    contig_length = contig_lengths[contig]
    if contig_length == 0:
        if "|" in contig:
            # ambig. assignment
            total = 0
            individual = []
            for ctg in contig.split("|"):
                ctg_len = contig_lengths[ctg]
                total += ctg_len
                individual.append(str(ctg_len))
            individual = "|".join(individual)
            contig_length = f"SUM:{total}|{individual}"
    # unfortunate consequence: numbers as strings
    return str(contig_length)


def main():

    args = parse_command_line()

    # record optional columns for output re-ordering
    opt_columns = []
    if args.fasta is not None:
        assert args.fasta.is_file()
        contig_lengths = read_sequence_lengths(args.fasta)
    else:
        opt_columns = opt_columns["contig_length_bp"]
        contig_lengths = None


    graph, node_lengths = build_graph(args.graph)
    # create tig mapping
    layout = read_layout_file(args.layout)

    unitigs = process_connected_components(
        graph, node_lengths, layout
    )

    hifi_cov = read_coverage_file(args.hifi_cov, "hifi")
    ont_cov = read_coverage_file(args.ont_cov, "ont")
    unitigs = unitigs.merge(hifi_cov, on="node", how="outer")
    unitigs["hifi_cov"].fillna(0, inplace=True)

    unitigs = unitigs.merge(ont_cov, on="node", how="outer")
    unitigs["ont_cov"].fillna(0, inplace=True)
    assert not pd.isnull(unitigs).all(axis=0).all()

    if contig_lengths is not None:
        # NB: use of collection.Counter() implied auto-default of 0
        assign_length = fnt.partial(assign_contig_length, contig_lengths)
        unitigs["contig_length_bp"] = unitigs["contig"].apply(assign_length)

    if args.rukki_paths is not None:
        paths, nodes_to_path_id = compute_rukki_path_ids(args.rukki_paths)
        unitigs["path_id"] = unitigs["unitig"].apply(
            lambda x: nodes_to_path_id.get(x, "no-path-member")
        )
        if args.path_ids is not None:
            args.path_ids.parent.mkdir(exist_ok=True, parents=True)
            paths.to_csv(args.path_ids, sep="\t", header=True, index=False)
    else:
        opt_columns.append("path_id")


    if args.rukki_final is not None:
        unitigs = process_final_node_assignment_table(unitigs, args.rukki_final)
    else:
        opt_columns.extend(["node_class", "class_color"])

    if args.rukki_colors is not None:
        unitigs = process_rukki_color_table(unitigs, args.rukki_colors)
    else:
        opt_columns.extend(
            ["rukki_color", "mat_marker", "pat_marker", "mat_div_pat"]
        )

    assert (unitigs["unitig"] == unitigs["node"]).all()
    unitigs.drop("unitig", axis=1, inplace=True)

    if "node_class" in unitigs.columns:
        unitigs.loc[unitigs["node_class"] == "maternal", "hap"] = "H1"
        unitigs.loc[unitigs["node_class"] == "paternal", "hap"] = "H2"

    order_columns = [
        "node", "hap", "mat_marker", "pat_marker",
        "rukki_color", "contig", "contig_length_bp",
        "node_class", "mat_div_pat",
        "hifi_cov", "ont_cov", "class_color", "path_id",
        "cc_id", "cc_size", "cc_length_hpc", "node_length_hpc",
        "layout_assigned", "is_used", "contig_pieces", "pieces_path",
        "total_gap_length", "default_gaps", "min_gaps",
        "ambiguous", "node_num"
    ]
    order_columns = [c for c in order_columns if c not in opt_columns]
    assert all([c in order_columns for c in unitigs.columns])

    unitigs = unitigs[order_columns]

    unitigs.sort_values(
        ["node_num", "hap", "node_length_hpc", "cc_id"],
        ascending=[True, False, False, False],
        inplace=True
    )

    args.table.parent.mkdir(exist_ok=True, parents=True)
    unitigs.to_csv(args.table, sep="\t", header=True, index=False)

    return


if __name__ == "__main__":
    main()
