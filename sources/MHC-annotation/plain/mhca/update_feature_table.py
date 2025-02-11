from argparse import ArgumentParser 
from collections import defaultdict
import mhca
import sys

def print_problem(prob):
    print(f"Problem found. {prob}")
    sys.exit(1)
    

class Gene:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
        if start < stop:
            self.strand = "+"
        else:
            self.strand = "-"
        self.sub_features = []
        self.short_name = ""
        self.locus_tag = ""
        self.type = ""
    def __str__(self):
        return f"{self.short_name}: {self.start}-{self.stop} ({self.strand}). Nr. of sub-features: {len(self.sub_features)}"
    def __lt__(self, other):
        return self.start < other.start

    def get_tbl_formated_coordinate_strings(self,c1,c2):
        f1 = (1 + (len(str(c1)) // 8)) * 8 
        f2 = (1 + (len(str(c2)) // 8)) * 8
        return f"{c1:<{f1}}", f"{c2:<{f2}}"

    def table_str(self, locus_tag=None):
        s1, s2 = self.get_tbl_formated_coordinate_strings(self.start,self.stop)
        outstr = f"{s1}{s2}gene\n"
        outstr += f"{' '*24}gene    {self.short_name}\n"
        if locus_tag: outstr += f"{' '*24}locus_tag       {locus_tag}\n"
        else: outstr += f"{' '*24}locus_tag       {self.locus_tag}\n"
        if self.type == "pseudogene":
            outstr += f"{' '*24}pseudogene      unknown\n"
        for subf in self.sub_features:
            for pidx, interval in enumerate(subf.intervals):
                s1, s2 = self.get_tbl_formated_coordinate_strings(interval[0], interval[1])
                if pidx == 0: outstr += f"{s1}{s2}{subf.type}\n"
                else: outstr += f"{s1}{s2.rstrip()}\n"
            outstr += f"{' '*24}product {self.short_name}\n"
            if subf.type == "CDS":
                outstr += f"{' '*24}transl_table    1\n"
            if subf.transcript_id != "": outstr += f"{' '*24}transcript_id   {subf.transcript_id}\n"
            else: outstr += f"{' '*24}transcript_id   {subf.ntranscript_id}\n"

        return outstr

class sub_feature:
    def __init__(self, start, stop, stype):
        self.start = start
        self.stop = stop
        if start < stop:
            self.strand = "+"
        else:
            self.strand = "-"
        self.intervals = [(start, stop)]
        self.transl_table = ""
        self.type = stype
        self.product = ""
        self.transcript_id = ""

    def __str__(self):
        return f"{self.short_name}: {self.start}-{self.stop} ({self.strand}). Nr. of features: {len(self.sub_features)}"
        
def parse_feature_table(feature_table_file, several=False):
    if several: all_features = defaultdict(set)
    else: features = set()
    with open(feature_table_file) as inf:
        seqid = ""
        current_feat = None
        for line_nr, line in enumerate(inf):
            #print(f"{line_nr}, length: {len(all_features['Feature gb|OK649231|'])}")
            if line.startswith(">"): 
                if not several: print_problem(f"Only one feature table expected but line with '>' found. {line}")
                if seqid: all_features[seqid].add(current_gene)
                seqid = line.lstrip(">").rstrip()
                current_feat = None
                continue
            else:
                if current_feat in {None, "end_sub", "g3"}:
                    if line[0] == " ": 
                        current_gene.type = "pseudogene"
                        if several: all_features[seqid].add(current_gene)
                        else: features.add(current_gene)
                        current_feat = None
                        continue
                    if len(line.split()) != 3: print_problem(f"feat: {current_feat}, {line}")
                    start, stop, rfeat = line.rstrip().split()
                    if rfeat == "gene": 
                        if current_feat: 
                            if several: all_features[seqid].add(current_gene)
                            else: features.add(current_gene)

                        current_gene = Gene(int(start), int(stop))
                        current_gene.type = "gene"
                        current_feat = "g1"
                    else:
                        current_sub = sub_feature(int(start), int(stop), rfeat)
                        if rfeat == "mRNA":
                            current_feat = "m1"
                        elif rfeat == "CDS":
                            current_feat = "c1"
                            current_sub.transl_table = "1"
                        current_gene.sub_features.append(current_sub)

                elif current_feat == "g1":
                    if line[0] != " ": print_problem(f"feature: {current_feat}")
                    gene = line.rstrip().split()[1]
                    current_gene.short_name = gene
                    current_feat = "g2"
                elif current_feat == "g2":
                    locus_tag = line.rstrip().split()[1]
                    current_gene.locus_tag = locus_tag
                    current_feat = "g3"
                elif current_feat == "m1":
                    if line[0] == " ":
                        current_sub.product = line.rstrip().split()[1]
                        current_feat = "m2"
                    else:
                        start, stop = line.rstrip().split()
                        current_sub.intervals.append((int(start), int(stop)))
                elif current_feat == "c1":
                    if line[0] == " ":
                        current_sub.product = line.rstrip().split()[1]
                        current_feat = "c2"
                    else:
                        start, stop = line.rstrip().split()
                        current_sub.intervals.append((int(start), int(stop)))
                elif current_feat == "c2":
                    current_feat = "c3"
                elif current_feat in {"c3","m2"}:
                    current_sub.transcript_id = line.rstrip().split()[1]
                    current_feat = "end_sub"
                        
                else:
                    print_problem(f"feat: {current_feat}")
        if several: 
            all_features[seqid].add(current_gene)
            return all_features
        else: 
            features.add(current_gene)
            return features

def main(args):
    """ Main function extract full fasta from refseq file"""

    allf = parse_feature_table(args.old_feature_table, several=True)
    #for sid, gs in allf.items(): 
        #for g in gs: print(g)
        #print(len(gs))
    #print(allf.keys())
    old_features = allf[args.record_id]
    old_feat_d = {}
    for feat in old_features:
        old_feat_d[feat.short_name] = feat
    new_features = parse_feature_table(args.update_feature_table, several=False)
    new_feat_d = {}
    for feat in new_features:
        new_feat_d[feat.short_name] = feat
    #for f in feats:
    #    print(f)
    print(f"Number of genes in old table: {len(old_feat_d)}")
    print(f"Number of genes in new table: {len(new_feat_d)}")
    print(f"Number of genes overlapping: {len(old_feat_d.keys() & new_feat_d.keys())}")
    #for feat in old_features-feats:
    #    print(feat)
        
    #print(f"Number of genes to be added: {len(feats-old_features)}")
    nrs = []
    for name, gene in old_feat_d.items():
        pref, nr = gene.locus_tag.split("_")
        nrs.append(int(nr))
    max_nr = max(nrs)
    cnr = max_nr
    for gname in new_feat_d.keys() - old_feat_d.keys():
        cnr += 1
        new_feat_d[gname].locus_tag = pref + "_" + f"{cnr:04d}"
        for sfeat in new_feat_d[gname].sub_features:
            sub_prefix, nr = sfeat.transcript_id.split("_")
            _, sub_nr = nr.split(".")
            sfeat.transcript_id = sub_prefix + f"_{cnr:04d}.{sub_nr}"
       
        
    with open(args.new_feature_table, "w") as outf:
        outf.write(f">{args.record_id}\n")
        for gname in sorted(new_feat_d, key = lambda x: new_feat_d[x].start):
            if gname in old_feat_d: 
                outf.write(old_feat_d[gname].table_str())
            else:
                outf.write(new_feat_d[gname].table_str())
    """
    for feat in old_features:
        if feat.short_name == "GPX5": 
            print(feat)
            print(feat.locus_tag)
    #print(feats)
    """

            
                
            
            
