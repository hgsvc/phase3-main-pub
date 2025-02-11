#!/usr/bin/env python3

# MELT-RISC: identify and classify mobile element insertions (MEIs) in PAV calls.

import argparse
import csv
import gzip
import hashlib
import os
import re
import sys
from ncls import NCLS

# ------------------------------------------------------
# Globals
# ------------------------------------------------------
VERSION = '1.4.1'
FASTA_SUFFIX_RE = r'\.(fa.gz|fasta.gz)$'
FASTA_FILE_RE = r'^(.*)' + FASTA_SUFFIX_RE
DEBUG = False

# increase CSV max field size
csv.field_size_limit(256 * 1024 * 1024)

MIN_TSD_LEN = 1
MAX_TSD_LEN = 40

# perl -e 'while (<>) { chomp; if (!(/^>/)) { $l += length($_);} } print "length=$l\n"' <SVA_A.fa 
ME_LENGTHS = {
    'ALU': 281,
    'SVA': 1316,    # MELT SVA reference
    'SVA_A': 1387,
    'SVA_F': 1375,
    'LINE1': 6019
}

# CSV output
CSV_HEADERS = ['samples', 'chrom', 'pos', 'strand', 'ME', '%ME', '%id', '%cov', 'insertion_seq', 'left_flank_seq', 'right_flank_seq', 'TSD_seq', 'polyX_coords', 'ME_coords', 'insertion_coords', 'match_string',
               'ME_family', 'ME_subfamily', 'ME_start', 'ME_stop', 'ME_num_diag_matches', 'ME_num_diffs', 'ME_diffs', 'overlapping_annots', 'genotype', 'hap1_region', 'hap2_region']
CSV_FLANKING_SEQ_BP = 30

# reverse complement
NA_CHARS = "actgnACTGN"
NA_COMP = "tgacnTGACN"
NA_RE = re.compile('^[' + NA_CHARS + ']*$')
NA_TRANS = str.maketrans(NA_CHARS, NA_COMP)

# ------------------------------------------------------
# logging
# ------------------------------------------------------
def fatal(msg):
    sys.stderr.write("FATAL - " + msg + "\n")
    sys.exit(1)

def info(msg):
    sys.stderr.write("INFO - " + msg + "\n")

def debug(msg):
    if DEBUG:
        sys.stderr.write("DEBUG - " + msg + "\n")

def warn(msg):
    sys.stderr.write("WARN - " + msg + "\n")
    
# ------------------------------------------------------
# read_single_fasta_file
# ------------------------------------------------------
def read_single_fasta_file(fpath):
    defline = None
    seq = ''
    with gzip.open(fpath, 'rt') as fh:
        for line in fh:
            # defline
            m = re.match(r'^>(.*)$', line)
            if m:
                if defline is not None:
                    fatal("multiple deflines found in single sequence FASTA file")
                defline = m.group(1)
            else:
                rsl = line.rstrip()
                seq = seq + rsl
                rsl = rsl.upper()
                rsl = re.sub(r'[\s\-]+', '', rsl)

    seqlen = len(seq)
    info("read sequence " + defline + " of length " + str(seqlen) + " from " + fpath)
    return { 'defline': defline, 'seq': seq, 'len': seqlen }

# ------------------------------------------------------
# read_fasta_dir
# ------------------------------------------------------
def read_fasta_dir(dpath, seqid, skip_seqids):
    files = {}

    for file in os.listdir(dpath):
        m = re.match(FASTA_FILE_RE, file)
        if m:
            # read specified reference sequence only
            if ((seqid is not None) and (m.group(1) != seqid)) or (m.group(1) in skip_seqids):
                continue
            file_id = m.group(1)
            fpath = os.path.join(dpath, file)
            ff = read_single_fasta_file(fpath)
            if file_id in files:
                fatal("duplicate FASTA file id " + file_id)
            files[file_id] = ff
    return files
        
# ------------------------------------------------------
# read_vcf_contigs
# ------------------------------------------------------
def read_vcf_contigs(vpath):
    contigs_l = []
    contigs_d = {}
    lnum = 0
    with gzip.open(vpath, 'rt') as fh:
        for line in fh:
            lnum += 1
            m = re.match(r'^\#\#contig=\<ID=([^,]+),length=(\d+),md5=([a-z0-9]{32})\>$', line)
            if m:
                c_id = m.group(1)
                c_len = int(m.group(2))
                c_md5 = m.group(3)
                debug("file id=" + c_id + " len=" + str(c_len) + " md5=" + c_md5)
                contig = { 'id': c_id, 'len': c_len, 'md5': c_md5 }
                contigs_l.append(contig)
                if id in contigs_d:
                    fatal("duplicate contig id " + c_id + " at line " + str(lnum))
                contigs_d[c_id] = contig

    n_contigs = len(contigs_l)
    msg ="read " + str(n_contigs) + " contig(s) from " + vpath 
    info(msg)
    return { 'contigs_l': contigs_l, 'contigs_d': contigs_d }

# ------------------------------------------------------
# read_vcf_insertions
# ------------------------------------------------------
def read_vcf_insertions(vpath, ref_seqs, seqid, skip_seqids):
    insertions = []
    lnum = 0

    with gzip.open(vpath, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            lnum += 1
            if re.match(r'^#', row[0]):
                continue
            
            (chrom, pos, vcf_id, ref, alt, qual, filt, info, fmt, gt) = row
            pos = int(pos)

            if ((seqid is not None) and (chrom != seqid)) or (chrom in skip_seqids):
                continue
            
            inf_d = {}
            for inf in info.split(';'):
                (k, v) = inf.split('=')
                if k in inf_d:
                    fatal("key " + k + " already seen in " + info)
                inf_d[k] = v

            if inf_d['SVTYPE'] != 'INS':
                continue

            # sanity check - reference sequence is a single base
            if len(ref) != 1:
                fatal("len(REF) != 1")

            # sanity check - reference base equals first base of alt sequence
            if alt[0] != ref:
                fatal("ALT[0] (" + alt[0] + " != REF (" + ref + ")")

            # sanity check - check reference base against reference sequence
            refseq = ref_seqs[chrom]
            ref_base = refseq['seq'][pos-1]
            if ref_base.upper() != ref:
                fatal("reference file base at " + str(pos) + " (" + ref_base  + ") != REF (" + ref + ")")

            insertion = {
                'chrom': chrom,
                'pos': pos,
                'vcf_id': vcf_id,
                'ref': ref,
                'alt': alt,
                'len': len(alt) - 1,
                'ins': alt[1:],
                'qual': qual,
                'filt': filt,
                'info': info,
                'fmt': fmt,
                'gt': gt,
            }

            # assign hap regions to hap1/hap2 based on genotype - for homozygous hap1 should come first
            gtypes = gt.split('|')
            hap_nums = []
            hap_num = 1
            for gtype in gtypes:
                if gtype == '1':
                    hap_nums.append(hap_num)
                hap_num += 1
                    
            # hap1_region, hap2_region
            hap_nums_ind = 0
            tig_region = inf_d['TIG_REGION']
            regions = tig_region.split(',')
            for region in regions: 
               m = re.match(r'^([^:]+:\d+-\d+)$', region)
               if m:
                   hap_region = m.group(1)
                   insertion["hap" + str(hap_nums[hap_nums_ind]) + "_region"] = hap_region
                   hap_nums_ind += 1
               else:
                   fatal("unable to parse TIG_REGION " + region + " at line " + str(lnum))

            insertions.append(insertion)
            
    return insertions

# ------------------------------------------------------
# read_ucsc_rmsk()
# ------------------------------------------------------
def read_ucsc_rmsk(ucsc_rmsk):
    annots_by_chrom = {}
    nr = 0
    starts = []
    ends = []
    ids = []
    
    with gzip.open(ucsc_rmsk, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            (bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id) = row
            chrom = genoName
            chromStart = genoStart
            chromEnd = genoEnd
            if chrom not in annots_by_chrom:
                annots_by_chrom[chrom] = { 'starts': [], 'ends': [], 'n': 0, 'rows': [], 'ncl': None }
            annots_by_chrom[chrom]['starts'].append(int(chromStart))
            annots_by_chrom[chrom]['ends'].append(int(chromEnd))
            annots_by_chrom[chrom]['n'] += 1
            annots_by_chrom[chrom]['rows'].append(row)
            nr += 1

    info("read " + str(nr) + " repeat rows from " + ucsc_rmsk)

    # create NCLs
    for chrom in annots_by_chrom:
        bc = annots_by_chrom[chrom]
        bc['ncl'] = NCLS(bc['starts'], bc['ends'], range(0, bc['n']))

    return annots_by_chrom
        
def get_overlapping_annotation(annots_by_chrom, chrom, start, end):
    # not every contig has repeatmasker annotations in UCSC db
    if chrom not in annots_by_chrom:
        return []
    by_chrom = annots_by_chrom[chrom]
    annots = []
    for j in by_chrom['ncl'].find_overlap(start, end):
        annots.append(by_chrom['rows'][j[2]])
    return annots

# ------------------------------------------------------
# read_water()
# ------------------------------------------------------
def read_water(wfile):
    n_lines = 0
    indel_lengths = {}
    alignment = None
    # alignments indexed by insertion seq name e.g., chr1-191377-INS-2012
    alignments = {}

    def process_alignment(al):
        vn = al['variant_name']
        if vn in alignments:
            fatal("duplicate variant name " + vn)
        alignments[vn] = alignment

        def update_coords(coords, mline):
            m = re.match(r'^\S+\s+(\d+) \S+\s+(\d+)\s*$', mline)
            if not m:
                fatal("unable to parse match line '" + mline + "'")
            c1 = int(m.group(1))
            c2 = int(m.group(2))
            (l, h) = (c1, c2) if c1 <= c2 else (c2, c1)
            if coords['x1'] is None or coords['x1'] > l:
                coords['x1'] = l
            if coords['x2'] is None or coords['x2'] < h:
                coords['x2'] = h
        
        # parse min/max alignment coords from match lines
        me_coords = { 'x1': None, 'x2': None }
        insertion_coords = { 'x1': None, 'x2': None }
        ctr = -1

        # create cigar/match string, using "^" for gap in ref seq, "v" for gap in insertion seq:
        #
        # input:
        #   ALU               55 --GGCGGATCACGAGGTCAGG---AGATCGAGACCATCCTGGCTAACACG     99
        #                          ||||           ||||   |||  |||.|...||..||     ||
        #   chr1-10829-IN     55 CCGGCG-----------CAGGCGCAGA--GAGGCGCGCCGCGC-----CG     86
        #
        # output:
        #   ^^||||vvvvvvvvvvv||||^^^|||vv|||.|...||..||vvvvv||
        # 
        al['match_str'] = ""
        match_chars = None
        seq_x1 = None
        seq_x2 = None
        
        def update_match_str(ct, ml):
            nonlocal match_chars
            nonlocal seq_x1
            nonlocal seq_x2
            
            if ct == 1 or ct == 3:
                m = re.match(r'^(\S+\s+\d+ )(\S+)(\s+\d+)\s*$', ml)
                if not m:
                    fatal("unable to parse match line 1 - " + ml)
                if ct == 1:
                    match_chars = list(re.sub(r'-', '^', m.group(2)))
                    seq_x1 = len(m.group(1))
                    seq_x2 = seq_x1 + len(match_chars)
                else:
                    seq = m.group(2)
                    for i in range(0, len(seq)):
                        if seq[i] == '-':
                            match_chars[i] = "v"
                    al['match_str'] = al['match_str'] + "".join(match_chars)
            elif ct == 2:
                matches = ml[seq_x1:seq_x2]
                for i in range(0, len(matches)):
                    if matches[i] != ' ':
                        match_chars[i] = matches[i]
        
        for ml in al['match_lines']:
            ctr = (ctr + 1) % 4
            if ctr == 1:
                if re.match(r'^\s*$', ml):
                    continue
                update_coords(me_coords, ml)
            elif ctr == 3:
                update_coords(insertion_coords, ml)
            update_match_str(ctr, ml)
        
        al['ME_coords'] = me_coords
        al['insertion_coords'] = insertion_coords
        # save memory
#        al['match_lines'] = None
        
    # reading detailed alignment match lines
    reading_match_lines = False
        
    with gzip.open(wfile, 'rt') as fh:
        for line in fh:
            n_lines += 1
            if re.match(r'^\# Aligned_sequences: 2.*$', line):
                if alignment is not None:
                    process_alignment(alignment)
                alignment = { 'match_lines': [] }
            else:
                m = re.match(r'^# 1: (\S+).*', line)
                if m:
                    alignment['ME'] = m.group(1)
                    
                m = re.match(r'^# 2: (\S+-(INS|DEL)-(\d+))', line)
                if m:
                    alignment['variant_name'] = m.group(1)
                    alignment['variant_type'] = m.group(2)
                    alignment['variant_len'] = int(m.group(3))
                    if alignment['variant_len'] < 50:
                        print("lnum=" + str(n_lines))

                m = re.match(r'^# Identity:\s+(\d+)\/(\d+) \((\s*[\d\.]+)%\)', line)
                if m:
                    alignment['matches'] = int(m.group(1))
                    alignment['length'] = int(m.group(2))
                    alignment['pct_id'] = float(m.group(3))

                m = re.match(r'^# Gaps:\s+(\d+)\/(\d+) \((\s*[\d\.]+)%\)', line)
                if m:
                    alignment['gaps'] = int(m.group(1))

                m = re.match(r'^# Score: ([\d\.]+).*$', line)
                if m:
                    alignment['score'] = float(m.group(1))

                m = re.match(r'^#=======================================\s*$', line)
                if alignment is not None and 'ME' in alignment and m:
                    reading_match_lines = True if not reading_match_lines else False
                elif re.match(r'^#\-+\s*$', line):
                    pass
                elif reading_match_lines:
                    alignment['match_lines'].append(line)
                    
    process_alignment(alignment)
    return alignments

# ------------------------------------------------------
# check_insertion_for_tsd
# ------------------------------------------------------
def check_insertion_for_tsd(ins, ref_seqs):
    ref_seq = ref_seqs[ins['chrom']]
    ref_seq_len = ref_seq['len']
    ins_seq = ins['ins']
    ins_seq_len = len(ins_seq)

    # find longest TSD _after_ the insertion

    # forward strand:
    # [<TSD>........polyA]<TSD>
    #
    # reverse strand:
    # [<TSD>polyT........]<TSD>
    #
    # [] = inserted sequence

    # note that the following situation is possible too, but it depends on the insertion caller
    # where the 5' TSD appears wrt to the called insertion. PAV always seems to place the
    # 5' TSD inside the insertion
    #
    # <TSD>][........polyA<TSD>]
    
    ref_seq_pos = ins['pos']
    ins_seq_pos = 0
    tsd_after = ''

    while ref_seq['seq'][ref_seq_pos].upper() == ins_seq[ins_seq_pos].upper():
        tsd_after = tsd_after + ref_seq['seq'][ref_seq_pos]
        ref_seq_pos += 1
        ins_seq_pos += 1
        if ref_seq_pos >= ref_seq_len or ins_seq_pos >= ins_seq_len:
            break

    # find longest TSD _before_ the insertion:
    # <TSD>][........polyA<TSD>]
    ref_seq_pos = ins['pos'] - 1
    ins_seq_pos = ins_seq_len - 1
    tsd_before = ''
    while ref_seq['seq'][ref_seq_pos].upper() == ins_seq[ins_seq_pos].upper():
        tsd_before = ref_seq['seq'][ref_seq_pos] + tsd_before
        ref_seq_pos -= 1
        ins_seq_pos -= 1
        if ref_seq_pos < 0 or ins_seq_pos < 0:
            break

    ins['tsds'] = {
        'before': { 'tsd': tsd_before, 'len': len(tsd_before), 'x1': ins_seq_len - len(tsd_before), 'x2': ins_seq_len - 1},
        'after': { 'tsd': tsd_after, 'len': len(tsd_after), 'x1': 1, 'x2': len(tsd_after)},
    }

def pfeat(l, x1, x2, name):
    fl = x2 - x1 + 1
    str = "".ljust(x1-1) + name.center(fl, "-") + "".ljust(l - x2 + 1)
    return str

def pfeat2(l, f, name):
    x1 = f['x1']
    x2 = f['x2']
    fl = x2 - x1 + 1
    str = "".ljust(x1-1) + "".center(fl, "-") + "".ljust(l - x2 + 1) + "  <- " + name
    return str

# Take the union of a set of intervals
def union_intervals(ivs):
    # sort by start coord
    sorted_ivs = sorted(ivs, key = lambda x: x['x1'], reverse=False)
    # join overlapping intervals
    i = 0
    while i < len(sorted_ivs) - 1:
        i1 = sorted_ivs[i]
        i2 = sorted_ivs[i+1]
        # overlap, merge and check for overlap with next
        if i1['x2'] >= i2['x1']:
            i1['x2'] = i2['x2']
            sorted_ivs.pop(i+1)
            # no overlap, consider next pair
        else:
            i += 1
    return sorted_ivs

# Intersect a list of intervals with a single specified interval
def intersect_intervals(ivs, ii):
    res_ivs = []
    for i in ivs:
        new_x1 = None
        new_x2 = None
        # ii starts inside i
        if i['x1'] <= ii['x1'] <= i['x2']:
            new_x1 = ii['x1']
            new_x2 = ii['x2'] if ii['x2'] <= i['x2'] else i['x2']
        # ii ends inside i
        elif i['x1'] <= ii['x2'] <= i['x2']:
            new_x1 = i['x1']
            new_x2 = ii['x2']
        # ii contains i
        elif ii['x1'] <= i['x1'] and ii['x2'] >= i['x2']:
            new_x1 = i['x1']
            new_x2 = i['x2']
        
        if new_x1 is not None:
            res_ivs.append({'x1': new_x1, 'x2': new_x2})
        
    return res_ivs

def sum_interval_lengths(ivs):
    s = 0
    for i in ivs:
        s += i['x2'] - i['x1'] + 1
    return s

# convert match string from read_water to list of alignment spans
def parse_alignment_spans(align, strand):
    ins_x1 = align['insertion_coords']['x1']
    ins_x2 = align['insertion_coords']['x2']
    me_x1 = align['ME_coords']['x1']
    me_x2 = align['ME_coords']['x2']
    match_str = align['match_str']
    msl = len(match_str)
    spans = []

    debug("ins coords=" + str([ins_x1, ins_x2]) + " me coords=" + str([me_x1, me_x2]))
    
    # offsets from left side of the alignment - will add ins_x1, me_x1 to these
    ins_o1 = 0
    ins_o2 = 0
    me_o1 = 0
    me_o2 = 0
    n_id_bp = 0

    # add alignment span
    def add_span():
        nonlocal ins_o1, ins_o2, me_o1, me_o2, n_id_bp
        if (((ins_o2 - ins_o1) != 0) and ((me_o2 - me_o1) != 0)):
            # handle reverse strand matches
            ins_c = None
            me_c = None
            len1 = None
            
            if (strand == '+'):
                ins_c = [ins_x1 + ins_o1 - 1, ins_x1 + ins_o2 - 1]
                me_c = [me_x1 + me_o1 - 1, me_x1 + me_o2 - 1]
                len1 = ins_c[1] - ins_c[0] + 1
                if ins_c[0] < 0:
                    fatal("ins_c[0] < 0, ins_c=" + str(ins_c))
            else:
                ins_c = [ins_x2 - ins_o1, ins_x2 - ins_o2]
                me_c = [me_x1 + me_o1 - 1, me_x1 + me_o2 - 1]
                len1 = ins_c[0] - ins_c[1] + 1
                if ins_c[1] < 0:
                    fatal("ins_c[1] < 0, ins_c=" + str(ins_c))
                
            span = {
                'ins': ins_c, 
                'me': me_c,
                'pct_id': (n_id_bp / len1) * 100.0
            }

            # sanity check
            len2 = span['me'][1] - span['me'][0] + 1
            if len1 != len2:
                debug("span=" + str(span))
                fatal("len1/len2 mismatch in add_span, len1=" + str(len1) + " len2=" + str(len2))
                
            spans.append(span)
            n_id_bp = 0

    # match string contains only ^ (gap in insertion), v (gap in ME), |, and .
    for i in range(0, msl):
        if (match_str[i] == '^'):
            add_span()
            ins_o2 += 1
            ins_o1 = ins_o2
            me_o1 = me_o2
        elif (match_str[i] == 'v'):
            add_span()
            me_o2 += 1
            ins_o1 = ins_o2
            me_o1 = me_o2
        else:
            if (match_str[i] == '|'):
                n_id_bp += 1
            # start or extend span
            ins_o2 += 1
            me_o2 += 1

    add_span()
    return spans

def get_insertion_alignments(vcf_ins, aligns):
    ins_name = vcf_ins['chrom'] + '-' + str(vcf_ins['pos']+1) + '-INS-' + str(vcf_ins['len'])
    if ins_name not in aligns:
        ins_name = vcf_ins['chrom'] + '-' + str(vcf_ins['pos']) + '-INS-' + str(vcf_ins['len'])
        if ins_name not in aligns:
            fatal("no alignment found for " + ins_name)
    al = aligns[ins_name]
    return (ins_name, al)

def trim_to_range(i, min, max):
    if i['x1'] < min:
        i['x1'] = min
    if i['x2'] > max:
        i['x2'] = max
    return i
    
def check_insertion_for_ME_match(vcf_ins, aligns, min_pctid, min_pctcov, strand):
    (ins_name, al) = get_insertion_alignments(vcf_ins, aligns)
    me_match = None

    # compute percent identity without gaps
    pctid_nogaps = (al['matches'] / (al['length'] - al['gaps'])) * 100.0
    
    debug("checking " + ins_name + " for match")
    if (pctid_nogaps >= min_pctid):
        debug("checking " + ins_name + " for match, pctid_nogaps is in range ")
        ins_len = len(vcf_ins['ins'])
        me_x1 = al['insertion_coords']['x1']
        me_x2 = al['insertion_coords']['x2']
        tsd = vcf_ins['tsds']['after'];

        # use orientation of match to break ties
        polyX = vcf_ins['polyA']
        if vcf_ins['polyT']['score'] > polyX['score']:
            polyX = vcf_ins['polyT']
        elif vcf_ins['polyT']['score'] == polyX['score']:
            polyX = vcf_ins['polyA'] if strand == '+' else vcf_ins['polyT']
        vcf_ins['polyX'] = polyX
            
        # coordinates should be base-based, indexed from 1
        if DEBUG:
            print(vcf_ins['ins'])
            print(pfeat(ins_len, me_x1, me_x2, al['ME']))
            if tsd['len'] > 0:
                print(pfeat(ins_len, tsd['x1'], tsd['x2'], 'TSD'))
            if polyX['len'] > 0:
                print(pfeat(ins_len, polyX['x1'], polyX['x2'], 'polyX'))

        # compute percent of the insertion covered by the ME alignment after removing TSD + polyX
        #  make list of intervals to remove (e.g., TSD, polyX)
        #  trim intervals to insertion coords (e.g., for super-long TSDs)
        intervals = [trim_to_range({'x1': i['x1'], 'x2': i['x2']}, 0, vcf_ins['len']) for i in [tsd, polyX] if i['len'] > 0]
        debug("insertion = " + vcf_ins['ins'])
        debug("intervals = " + str(intervals))
        # take the union in case TSD, polyX overlap
        joined_intervals = union_intervals(intervals)
        joined_intervals_bp = sum_interval_lengths(joined_intervals)
        # number of bp left in insertion after removing TSD, polyX
        remaining_bp = vcf_ins['len'] - joined_intervals_bp
        debug("joined_intervals = " + str(joined_intervals))
        debug("joined_intervals_bp = " + str(joined_intervals_bp))

        # extract alignment spans from match_str
        spans = parse_alignment_spans(al, strand)
        debug(ins_name + " spans: " + str(spans) + "\n")

        # compute precise %identity and %coverage from alignment spans
        total_span_bp = 0
        total_aligned_span_bp = 0
        
        for span in spans:
            sx1 = span['ins'][0]
            sx2 = span['ins'][1]
            if (sx1 > sx2):
                tmp = sx2
                sx2 = sx1
                sx1 = tmp
            span_bp = sx2 - sx1
                
            # subtract joined_intervals from span after converting to 1-based base coordinates
            isecs = intersect_intervals(joined_intervals, { 'x1': sx1+1, 'x2': sx2 })
            debug("intersecting joined intervals with " + str([sx1+1,sx2]) + " = " + str(isecs))
            
            # update total_span_bp, total_aligned_span_bp
            subtract_bp = 0
            for isec in isecs:
                # intersected spans use 1-based base coordinates
                diff = isec['x2'] - isec['x1'] + 1
                subtract_bp += diff

            if subtract_bp > span_bp:
                fatal("span_bp=" + str(span_bp) + " subtract_bp=" + str(subtract_bp))
            new_span_bp = span_bp - subtract_bp
            total_span_bp = total_span_bp + new_span_bp
            # TODO - note that this is an estimate if subtract_bp > 0: it assumes the subtracted region has the same
            # average percent identity as the entire alignment
            total_aligned_span_bp += (new_span_bp * (span['pct_id']/100.0))
                
        if total_span_bp > remaining_bp:
            fatal("total_span_bp (" + str(total_span_bp) + ") > remaining_bp (" + str(remaining_bp) + ")")

        span_rem_ins_pctcov = (total_span_bp / remaining_bp) * 100.0 if remaining_bp > 0 else 0
        span_rem_ins_pctid = (total_aligned_span_bp / total_span_bp) * 100.0 if total_span_bp > 0 else 0

        #  intersect overlapping intervals with alignment span
        intersected_intervals = intersect_intervals(joined_intervals, al['insertion_coords'])
        intersected_intervals_bp = sum_interval_lengths(intersected_intervals)
        debug("intersected_intervals = " + str(intersected_intervals))
        debug("intersected_intervals_bp = " + str(intersected_intervals_bp))
        
        # alignment length in insertion and ME coords
        al_len = al['insertion_coords']['x2'] - al['insertion_coords']['x1'] + 1
        al_len_me = al['ME_coords']['x2'] - al['ME_coords']['x1'] + 1
        
        # bp in the alignment not covered by TSD/polyX
        remaining_alignment_bp = al_len - intersected_intervals_bp
        rem_ins_pctcov = (remaining_alignment_bp / remaining_bp) * 100.0 if remaining_bp > 0 else 0
        debug("remaining_bp = " + str(remaining_bp) + " remaining_alignment_bp = " + str(remaining_alignment_bp))

        # sanity checks
        if span_rem_ins_pctcov < 0 or span_rem_ins_pctcov > 100:
            fatal("span_rem_ins_pctcov=" + str(span_rem_ins_pctcov))

        if span_rem_ins_pctid < 0 or span_rem_ins_pctid > 100:
            fatal("span_rem_ins_pctid=" + str(span_rem_ins_pctid))
        
        # percent of the _entire_ insertion covered by the ME alignment
        ins_pctcov = (al_len / vcf_ins['len']) * 100.0
        me_pctcov =  (al_len_me / ME_LENGTHS[al['ME']]) * 100.0
        debug("checking " + ins_name + " for match, span_rem_ins_pctcov=" + str(span_rem_ins_pctcov) + " me_pctcov=" + str(me_pctcov))

        if span_rem_ins_pctcov >= min_pctcov and span_rem_ins_pctid >= min_pctid:
            debug("checking " + ins_name + " for match, ins_pctcov is good, setting match to nonempty")
            me_match = {
                'ME': al['ME'],
                'pctid': span_rem_ins_pctid,
                'rem_ins_pctcov': span_rem_ins_pctcov,
                'ins_pctcov': ins_pctcov,
                'me_pctcov': me_pctcov,
                'alignment': al,
                'strand': strand
            }
    return me_match

def check_insertion_for_5_prime_inverted_ME_match(vcf_ins, aligns, rev_aligns, min_pctid, min_pctcov):
    (fwd_name, fwd_al) = get_insertion_alignments(vcf_ins, aligns)
    (rev_name, rev_al) = get_insertion_alignments(vcf_ins, rev_aligns)

    debug("checking " + fwd_name + " for 5'-inverted match")
    
    # do rough check of total coverage - we already know neither strand alone is sufficient
    fwd_x1 = fwd_al['insertion_coords']['x1']
    fwd_x2 = fwd_al['insertion_coords']['x2']

    rev_x1 = fwd_al['insertion_coords']['x1']
    rev_x2 = fwd_al['insertion_coords']['x2']

    #  trim intervals to insertion coords (e.g., for super-long TSDs)
    intervals = []
    if 'tsds' in vcf_ins:
        if 'after' in vcf_ins['tsds']:
            intervals.append(vcf_ins['tsds']['after'])

    if 'polyX' in vcf_ins:
        intervals.append(vcf_ins['polyX'])

    intervals = [trim_to_range({'x1': i['x1'], 'x2': i['x2']}, 0, vcf_ins['len']) for i in intervals if i['len'] > 0]
    debug("insertion = " + vcf_ins['ins'])
    debug("intervals = " + str(intervals))
    # take the union in case TSD, polyX overlap
    joined_intervals = union_intervals(intervals)
    joined_intervals_bp = sum_interval_lengths(joined_intervals)
    # number of bp left in insertion after removing TSD, polyX

#    remaining_bp = vcf_ins['len'] - joined_intervals_bp
    remaining_bp = vcf_ins['len']

    debug("joined_intervals = " + str(joined_intervals))
    debug("joined_intervals_bp = " + str(joined_intervals_bp))

    # TODO - subtract joined intervals from fwd/rev matches
    combined = union_intervals([{'x1': fwd_x1, 'x2': fwd_x2}, {'x1': rev_x1, 'x2': rev_x2}])
    combined_bp = 0
    for i in combined:
        combined_bp += i['x2'] - i['x1'] + 1

    me_match = None
    pctid_nogaps = (fwd_al['matches'] / (fwd_al['length'] - fwd_al['gaps'])) * 100.0
    pctcov =  (combined_bp / remaining_bp) * 100.0

    fwd_pctid = (fwd_al['matches'] / (fwd_al['length'] - fwd_al['gaps'])) * 100.0
    fwd_pctid_str = "%.1f" % fwd_pctid
    rev_pctid = (rev_al['matches'] / (rev_al['length'] - rev_al['gaps'])) * 100.0
    rev_pctid_str = "%.1f" % rev_pctid
    fwd_bp = fwd_x2 - fwd_x1 + 1
    rev_bp = fwd_x2 - fwd_x1 + 1
    avg_pctid = (fwd_pctid + rev_pctid) / 2
    avg_pctid_str = "%.1f" % avg_pctid

    if pctcov >= min_pctcov and avg_pctid > min_pctid:
        debug("fwd_bp=" + str(fwd_bp) + " fwd_pctid=" + fwd_pctid_str + " rev_bp=" + str(rev_bp) + " rev_pctid=" + rev_pctid_str)
        debug("avg_pctid=" + avg_pctid_str + " combined pctcov=" + str(pctcov))

        me_fwd_match = {
            'ME': fwd_al['ME'],
            'pctid': fwd_pctid,
            'rem_ins_pctcov': pctcov,
            'ins_pctcov': pctcov,
            'me_pctcov': pctcov,
            'alignment': fwd_al,
            'strand': '+'
        }

        me_rev_match = {
            'ME': rev_al['ME'],
            'pctid': rev_pctid,
            'rem_ins_pctcov': pctcov,
            'ins_pctcov': pctcov,
            'me_pctcov': pctcov,
            'alignment': rev_al,
            'strand': '-'
        }

        return (me_fwd_match, me_rev_match)
    return None

# find polyA/polyT using exact match
#
#  seq - sequence to search for polyA/polyT
#  base - either 'A' or 'T'
#  start - starting position betwee 0 and len(seq) - 1
#  offset - 1 to search forwards from start, -1 to search backwards from start
#
def find_polyX_exact(seq, base, start, offset):
    uc_base = base.upper()
    end = start
    plen = 0
        
    while (end >= 0) and (end < len(seq)) and seq[end].upper() == uc_base:
        end += offset
        plen += 1
        
    if offset > 0:
        x1 = start + 1
        x2 = end
    else:
        x1 = end + 2
        x2 = start + 1

    return { 'len': plen, 'x1': x1, 'x2': x2 }

# find polyA/polyT using sliding window
#
#  seq - sequence to search for polyA/polyT
#  base - either 'A' or 'T'
#  start - starting position betwee 0 and len(seq) - 1
#  offset - 1 to search forwards from start, -1 to search backwards from start
#  wlen - window length >= 1
#  max_mismatches - max number of bp mismatch in any window of length wlen
#
def find_polyX_sliding_window(seq, base, start, offset, wlen, max_mismatches):
    uc_base = base.upper()
    orig_start = start
    end = start
    base_count = 0
    nonbase_count = 0
    polyX_seq = ''

    debug("find_polyX_sliding_window base=" + str(base) + " offset=" + str(offset) + " start=" + str(start))
    
    def add_or_remove_base(pos, add, offset):
        nonlocal base_count
        nonlocal nonbase_count
        nonlocal polyX_seq
        
        if seq[pos].upper() == uc_base:
            base_count += add
        else:
            nonbase_count += add

        if add == 1:
            if offset > 0:
                polyX_seq = polyX_seq + seq[pos]
            else:
                polyX_seq = seq[pos] + polyX_seq

    # end is out of range
    if end >= len(seq):
        if offset == -1:
            return {'len': 0, 'x1': end+2, 'x2': end+1, 'score': 0}
        else:
            return {'len': 0, 'x1': end+1, 'x2': end, 'score': 0}
    
    # start with window of size 1
    wsize = 1
    add_or_remove_base(end, 1, offset)
        
    while ((offset == 1) or (end > 0)) and ((offset == -1) or (end < len(seq) - 1)) and (nonbase_count <= max_mismatches):
        # add a base to the end of the window
        end += offset
        wsize += 1
        add_or_remove_base(end, 1, offset)
            
        # remove a base from the start of the window
        if wsize > wlen:
            add_or_remove_base(start, -1, offset)
            start += offset
            wsize -= 1

    debug("done, end = " + str(end))
    debug("before nonbase trim = " + str(polyX_seq))
            
    # trim nonbase characters, but only from the end (start is given and assumed to be correct?)
    nb_re = '[^' + base + ']+'
    rep = nb_re + '$' if offset > 0 else '^' + nb_re
    new_polyX_seq = re.sub(rep, '', polyX_seq, re.IGNORECASE)
    bp_trimmed = len(polyX_seq) - len(new_polyX_seq)

    debug(" after nonbase trim = " + str(new_polyX_seq))
    debug("bp_trimmed=" + str(bp_trimmed))
    
    if offset > 0:
        x1 = orig_start + 1
        x2 = end + 1 - bp_trimmed
    else:
        x1 = end + 1 + bp_trimmed
        x2 = orig_start + 1

    debug("final coords = " + str(x1) + " - " + str(x2))
    # sanity check
    clen = x2 - x1 + 1
    if clen != len(new_polyX_seq):
        fatal("find_polyX computed polyX len=" + str(clen) + " actual length=" + str(len(new_polyX_seq)))
    
    only_matches_seq = re.sub(nb_re, '', new_polyX_seq, re.IGNORECASE)
        
    return { 'len': x2 - x1 + 1, 'x1': x1, 'x2': x2, 'score': len(only_matches_seq) }

def check_insertion_for_polyA(vcf_ins, window_size_bp, max_mismatch_bp):
    ins = vcf_ins['ins']
    # 1. search backwards from end of insertion sequence for polyA
#    vcf_ins['polyA'] = find_polyX_exact(ins, 'A', len(ins) - 1, -1)
    vcf_ins['polyA'] = find_polyX_sliding_window(ins, 'A', len(ins) - 1, -1, window_size_bp, max_mismatch_bp)
        
    # 2. search forwards from start of insertion sequence for polyT
#    vcf_ins['polyT'] = find_polyX_exact(ins, 'T', vcf_ins['tsds']['after']['len'], 1)
    vcf_ins['polyT'] = find_polyX_sliding_window(ins, 'T', vcf_ins['tsds']['after']['len'], 1, window_size_bp, max_mismatch_bp)

# ------------------------------------------------------
# write_me_fasta_fwd
# ------------------------------------------------------
def revcomp_na(seq):
    if not re.match(NA_RE, seq):
        fatal("NA sequence contains illegal character(s): " + seq)
    revcomp = seq[::-1].translate(NA_TRANS)
    return revcomp
    
# write ME sequence to FASTA file in forward orientation
def write_me_fasta_fwd(ins, fasta_fh):
    seqid = ins['chrom'] + ':' +  str(ins['pos']) + ':' + ins['me_match']['strand']
    
    # NOTE: we're passing the entire insertion seq but could limit it to the region that matches ref ME
    # shouldn't be much difference if everything's working correctly, and passing the entire seq makes
    # it easier to detect discrepancies between the two Smith-Waterman alignments (pipeline and cALU/LINEu)
    seq = ins['alt'][1:]

    # revcomp sequence if needed
    if ins['me_match']['strand'] == '-':
        seq = revcomp_na(seq)

    # write entry to FASTA
    fasta_fh.write(">" + seqid + "\n")
    fasta_fh.write(seq + "\n")
    return seqid

# ------------------------------------------------------
# run_melt_calu_lineu
# ------------------------------------------------------
def run_melt_calu_lineu(MEIs_d, melt_jar, melt_exec, fasta_file_path):
    if not re.match(r'^(CALU|LINEU)$', melt_exec):
        fatal("Unsupported MELT executable - " + melt_exec)
    melt_cmd = "java -jar " + melt_jar + " " + melt_exec + " -f " + fasta_file_path 
    melt_fh = os.popen(melt_cmd)
    reading = False
    
    for line in melt_fh:
        if re.match(r'^Performing MELT analysis.*$', line):
            reading = True
        elif re.match('^End time.*$', line):
            reading = False
        elif reading == True:
            m = re.match(r'^([^:]+:\d+:[\-\+])\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\S+|No Differences)\t(\S+)$', line)
            if not m:
                fatal("unable to parse MELT output: " + line)
            fasta_id = m.group(1)
            family = m.group(2)
            subfamily = m.group(3)
            start = int(m.group(4))
            stop =int(m.group(5))
            diag_matches =int(m.group(6))
            total_diffs =int(m.group(7))
            diffs = m.group(8)
            seq = m.group(9)

            mei = MEIs_d[fasta_id]
            if mei['fasta_id'] != fasta_id:
                fatal("seqid mismatch for " + mei['fasta_id'] + " / " + fasta_id)

            mei['ME_family'] = family
            mei['ME_subfamily'] = subfamily
            mei['ME_start'] = str(start)
            mei['ME_stop'] = str(stop)
            mei['ME_num_diag_matches'] = str(diag_matches)
            mei['ME_num_diffs'] = str(total_diffs)
            mei['ME_diffs'] = diffs
    
# ------------------------------------------------------
# add_insertion()
# ------------------------------------------------------
def ins_position_str(ins):
    return (ins['chrom'] + ":" + str(ins['pos'])).ljust(16)

def add_insertion(ins, ref_seqs, alu_fasta_fh, line_fasta_fh):
    ref_seq = ref_seqs[ins['chrom']]
    ref_seq_len = ref_seq['len']
    ref_seq_pos = ins['pos']
    l_context_bp = 10
    r_context_bp = 30

    def get_l_context_bp(cbp):
        bp_before = cbp if ref_seq_pos > cbp else ref_seq_pos
        sb_from = ref_seq_pos - bp_before
        return ref_seq['seq'][sb_from:ref_seq_pos]

    def get_r_context_bp(cbp):
        bp_after = cbp if (ref_seq_pos + cbp < ref_seq_len) else ref_seq_len - ref_seq_pos
        sa_to = ref_seq_pos + bp_after
        return ref_seq['seq'][ref_seq_pos:sa_to]

    # display region around insertion point (i.e., the point after the REF base)
    seq_before = get_l_context_bp(l_context_bp)
    seq_after = get_r_context_bp(r_context_bp)
    
    # last base of seq_before should be REF
    if seq_before[-1].upper() != ins['ref']:
        fatal("seq_before[-1] (" + seq_before[-1] + ") != REF (" + ins['ref'] + ")")

    tsd_str = ins['tsds']['after']['tsd'] if ins['tsds']['after']['tsd'] is not None else '-'
    pos_str = ins_position_str(ins)
    # sequence to display inside the insertion
    l_bp = 45
    r_bp = 45
    middle_bp = 14
    remaining_bp = ins['len'] - l_bp - r_bp
    ins_str = ins['alt'][1:l_bp+1] + "..." + ("+" + str(remaining_bp) + "bp").center(middle_bp-6,".") + "..." + ins['alt'][-r_bp:]

    pctid_str = (('%.1f' % ins['me_match']['pctid']) + "%").rjust(6)
    # percent of insert covered by alignment
    rem_ins_pctcov_str = (("%.1f" % ins['me_match']['rem_ins_pctcov']) + "%").rjust(6)
    # percent of reference ME sequence covered by alignment
    me_pctcov_str = (("%.1f" % ins['me_match']['me_pctcov']) + "%").rjust(6)
    me_str = ins['me_match']['ME'].ljust(5) + "|" + ins['me_match']['strand'].ljust(3) + "|" + me_pctcov_str + "|" + pctid_str + "|" + rem_ins_pctcov_str
    print(pos_str + "|" + me_str + "| " + seq_before + " [" + ins_str + "] " + seq_after)

    # print TSD, polyA position
    pos_str = re.sub(r'\S', ' ', pos_str)
    me_str = re.sub(r'\S', ' ', me_str)
    seq_before = re.sub(r'\S', ' ', seq_before)
    seq_after = ""
    isl = len(ins_str)
    ins_str_left = "".center(l_bp)
    ins_str_right = "".center(r_bp)
    ins_str_middle = "".center(middle_bp)
    
    # TSD
    tl = ins['tsds']['after']['len']
    if tl > 0:
        rep = '^' + 'TSD'.center(tl-2, '^') + '^'
        rep = rep[0:tl]
        ins_str_left = rep.ljust(l_bp)[0:l_bp]
        seq_after = rep

    # polyA
    pa = ins['polyA']
    pa_len = pa['len']
    if pa_len >= 7:
        ins_str_right = ('<' + "polyA".center(pa_len-2, "-") + '>').rjust(r_bp)[-r_bp:]
        
    # polyT
    pt = ins['polyT']
    pt_len = pt['len']
    if pt_len >= 7:
        ins_str_left = ins_str_left[0:pt['x1']-1] + '<' + "polyT".center(pt_len-2, "-") + '>' + ins_str_left[pt['x2']:]
        ins_str_left = ins_str_left[0:l_bp]
        
    print(pos_str + " " + me_str + "  " + seq_before + "  " + ins_str_left + ins_str_middle + ins_str_right + "  " + seq_after)

    # print ME alignment
    ins_coords = ins['me_match']['alignment']['insertion_coords']
    x1 = ins_coords['x1']
    x2 = ins_coords['x2']
    ins_ld = "<" if ins['me_match']['strand'] == '-' else '['
    ins_rd = "]" if ins['me_match']['strand'] == '-' else '>'
    ins_len = len(ins['ins'])
    ins_left = "".center(x1-1) + ins_ld + ins['me_match']['ME'].ljust(x2 - x1 + 1, "-") + ins_rd if x1 < l_bp else "".center(l_bp)
    ins_right = ins_ld + ins['me_match']['ME'].rjust(x2 - x1 + 1, "-") + ins_rd + "".center(ins_len - x2) if x2 > (ins_len - r_bp) else "".center(r_bp)
    ins_str = ins_left[0:l_bp] + "".center(14) + ins_right[-r_bp:]
    print(pos_str + " " + me_str + "  " + seq_before + "  " + ins_str)
    print()

    def coords2str(c):
        return str(c['x1']) + '-' + str(c['x2'])

    # ALU/LINE1 FASTA output
    fasta_id = None
    if (alu_fasta_fh is not None) and (ins['me_match']['ME'] == 'ALU'):
        fasta_id = write_me_fasta_fwd(ins, alu_fasta_fh)

    if (line_fasta_fh is not None) and (ins['me_match']['ME'] == 'LINE1'):
        fasta_id = write_me_fasta_fwd(ins, line_fasta_fh)
        
    l_flank = get_l_context_bp(CSV_FLANKING_SEQ_BP)
    r_flank = get_r_context_bp(CSV_FLANKING_SEQ_BP)
    ins_d = {
        'fasta_id': fasta_id,
        'chrom': ins['chrom'],
        'pos': str(ins['pos']),
        'strand': ins['me_match']['strand'],
        'ME': ins['me_match']['ME'],
        '%ME': me_pctcov_str.strip(),
        '%id': pctid_str.strip(),
        '%cov': rem_ins_pctcov_str.strip(),
        'insertion_seq': ins['alt'][1:],
        'left_flank_seq': l_flank,
        'right_flank_seq': r_flank,
        'TSD_seq': ins['tsds']['after']['tsd'],
        'polyX_coords': coords2str(ins['polyX']),
        'ME_coords': coords2str(ins['me_match']['alignment']['ME_coords']),
        'insertion_coords': coords2str(ins['me_match']['alignment']['insertion_coords']),
        'alignment': ins['me_match']['alignment']['match_str'],
        'genotype': ins['gt'],
        'hap1_region': ins['hap1_region'] if 'hap1_region' in ins else '',
        'hap2_region': ins['hap2_region'] if 'hap2_region' in ins else ''
    }
    
    return ins_d
    
# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():
    # input
    parser = argparse.ArgumentParser(description='Identify MEIs in a PAV-generated VCF file of sequence variants.')
    parser.add_argument('--vcf', required=True, help='Path to VCF file containing contigs to check.')
    parser.add_argument('--sample', required=True, help='Sample ID.')
    parser.add_argument('--fasta_dir', required=True, help='Path to directory that contains FASTA reference files.')
    parser.add_argument('--alu_water', required=True, help='Path to ALU water alignment output file.')
    parser.add_argument('--alu_water_rev', required=True, help='Path to ALU water reverse strand alignment output file.')
    parser.add_argument('--sva_water', required=True, help='Path to SVA water alignment output file.')
    parser.add_argument('--sva_water_rev', required=True, help='Path to SVA water reverse strand alignment output file.')
    parser.add_argument('--line_water', required=True, help='Path to LINE1 water alignment output file.')
    parser.add_argument('--line_water_rev', required=True, help='Path to LINE1 water reverse strand alignment output file.')
    parser.add_argument('--min_seqlen', required=False, type=int, default=100, help='Minimum insertion sequence length.')
    parser.add_argument('--max_seqlen', required=False, type=int, default=50000, help='Maximum insertion sequence length.')
    parser.add_argument('--min_pctid', required=False, type=int, default=90, help='Minimum average percent identity of aligned regions.')
    parser.add_argument('--min_pctcov', required=False, type=int, default=85, help='Minimum percent coverage of insertion minus TSD and polyX by aligned regions.')
    parser.add_argument('--polyx_window_bp', required=False, type=int, default=4, help='Sliding window size for polyA/polyT detection.')
    parser.add_argument('--polyx_max_mismatch_bp', required=False, type=int, default=1, help='Maximum number of mismatches in sliding window for polyA/polyT detection.')
    parser.add_argument('--seqid', required=False, help='Optional sequence id: process only insertions on this reference sequence.')
    parser.add_argument('--skip_seqids', required=False, help='Optional comma-delimited list of sequence ids to skip.')
    parser.add_argument('--csv_output', required=False, help='Optional path to CSV format output file.')
    parser.add_argument('--alu_fasta', required=True, help='Path to FASTA output file of forward-strand ALU sequences.')
    parser.add_argument('--line_fasta', required=True, help='Path to FASTA output file of forward-strand LINE1 sequences.')
    parser.add_argument('--melt_jar', required=True, help='Path to MELT JAR file for running CALU and LINEU subfamily analysis.')
    parser.add_argument('--ucsc_rmsk', required=False, help='Optional path to UCSC rmsk.txt.gz file containing repeat annotations.')
    args = parser.parse_args()

    skip_seqids = {}
    if args.skip_seqids is not None and args.skip_seqids != '':
        for seqid in args.skip_seqids.split(','):
            info("skipping insertions on sequence " + seqid)
            skip_seqids[seqid] = True

    # echo arguments
    info("VERSION " + VERSION)
    for arg in vars(args):
        info(arg + "=" + str(getattr(args, arg)))

    # index UCSC annotation, if present
    annots_by_chrom = None
    if args.ucsc_rmsk:
        annots_by_chrom = read_ucsc_rmsk(args.ucsc_rmsk)

    # read reference FASTA files
    fasta_files = read_fasta_dir(args.fasta_dir, args.seqid, skip_seqids)

    # read VCF contigs
    vcf_contigs = read_vcf_contigs(args.vcf)

    # TODO - incorporate MD5 checks from checkVCFReferenceSeqs.py
    
    # read VCF insertions
    vcf_insertions = read_vcf_insertions(args.vcf, fasta_files, args.seqid, skip_seqids)
    info("read " + str(len(vcf_insertions)) + " insertions from " + args.vcf)

    # filter by length
    vcf_insertions = [i for i in vcf_insertions if i['len'] >= args.min_seqlen and i['len'] <= args.max_seqlen]
    info("read " + str(len(vcf_insertions)) + " insertions with length >= " + str(args.min_seqlen) + " and length <= " + str(args.max_seqlen))
        
    # read alignment files
    alu_aligns = read_water(args.alu_water)
    info("read " + str(len(alu_aligns)) + " ALU alignment(s) from " + args.alu_water)
    alu_rev_aligns = read_water(args.alu_water_rev)
    info("read " + str(len(alu_rev_aligns)) + " reverse strand ALU alignment(s) from " + args.alu_water_rev)

    sva_aligns = read_water(args.sva_water)
    info("read " + str(len(sva_aligns)) + " SVA alignment(s) from " + args.sva_water)
    sva_rev_aligns = read_water(args.sva_water_rev)
    info("read " + str(len(sva_rev_aligns)) + " reverse strand SVA alignment(s) from " + args.sva_water_rev)

    line_aligns = read_water(args.line_water)
    info("read " + str(len(line_aligns)) + " LINE1 alignment(s) from " + args.line_water)
    line_rev_aligns = read_water(args.line_water_rev)
    info("read " + str(len(line_rev_aligns)) + " reverse strand LINE1 alignment(s) from " + args.line_water_rev)

    # FASTA output for CALU
    alu_fasta_fh = open(args.alu_fasta, "w")
    if alu_fasta_fh is None:
        fatal("failed to open/write Alu FASTA file " + args.alu_fasta)
    
    # FASTA output for LINEU
    line_fasta_fh = open(args.line_fasta, "w")
    if line_fasta_fh is None:
        fatal("failed to open/write LINE1 FASTA file " + args.line_fasta)

    # MEIs to print/report
    MEIs = []
    MEIs_d = {}
    
    # stdout summary output
    print("Location        |ME   |+/-|%ME   |%id   |%cov  | insertion")
        
    n_tsd_before = 0
    n_tsd_after = 0
    n_me_match = 0
    for vcf_ins in vcf_insertions:
        # check for presence of target site duplication
        check_insertion_for_tsd(vcf_ins, fasta_files)
        tsds = vcf_ins['tsds']

        if tsds['before']['len'] >= MIN_TSD_LEN and tsds['before']['len'] <= MAX_TSD_LEN:
            n_tsd_before += 1
        if tsds['after']['len'] >= MIN_TSD_LEN and tsds['after']['len'] <= MAX_TSD_LEN:
            n_tsd_after += 1

        # check for polyA/polyT
        check_insertion_for_polyA(vcf_ins, args.polyx_window_bp, args.polyx_max_mismatch_bp)

        if vcf_ins['polyA']['score'] > 3 and vcf_ins['polyT']['score'] > 3 and vcf_ins['polyA']['score'] == vcf_ins['polyT']['score']:
            warn(ins_position_str(vcf_ins) + " polyA and polyT both have score " + str(vcf_ins['polyA']['score']))
        
        # check for qualifying match with mobile element
        alu_match = check_insertion_for_ME_match(vcf_ins, alu_aligns, args.min_pctid, args.min_pctcov, '+')
        alu_rev_match = check_insertion_for_ME_match(vcf_ins, alu_rev_aligns, args.min_pctid, args.min_pctcov, '-')
        sva_match = check_insertion_for_ME_match(vcf_ins, sva_aligns, args.min_pctid, args.min_pctcov, '+')
        sva_rev_match = check_insertion_for_ME_match(vcf_ins, sva_rev_aligns, args.min_pctid, args.min_pctcov, '-')
        line_match = check_insertion_for_ME_match(vcf_ins, line_aligns, args.min_pctid, args.min_pctcov, '+')
        line_rev_match = check_insertion_for_ME_match(vcf_ins, line_rev_aligns, args.min_pctid, args.min_pctcov, '-')
        matches = [m for m in [alu_match, alu_rev_match, sva_match, sva_rev_match, line_match, line_rev_match] if m is not None]
        n_matches = len(matches)

        # sort matches and pick the highest-scoring
        me_match = None
        inv_matches = None

        if n_matches > 0:
            sorted_matches = sorted(matches, key = lambda x: x['alignment']['score'], reverse=True)
            me_match = sorted_matches[0]
            # kick the can down the road
            if n_matches > 1 and sorted_matches[0]['alignment']['score'] == sorted_matches[1]['alignment']['score']:
                warn("multiple qualifying matches with the same score: " + str(sorted_matches))

        # check for possible 5' inversions - do any matches qualify if both strands are combined?
        elif False:
            # TODO - may not be an issue for Alus (or SVA) due to >= 500bp requirement
            alu_matches = check_insertion_for_5_prime_inverted_ME_match(vcf_ins, alu_aligns, alu_rev_aligns, args.min_pctid, args.min_pctcov)
            sva_matches = check_insertion_for_5_prime_inverted_ME_match(vcf_ins, sva_aligns, sva_rev_aligns, args.min_pctid, args.min_pctcov)
            line1_matches = check_insertion_for_5_prime_inverted_ME_match(vcf_ins, line_aligns, line_rev_aligns, args.min_pctid, args.min_pctcov)

            if line1_matches is not None:
                inv_matches = line1_matches
            elif alu_matches is not None:
                inv_matches = alu_matches
            elif sva_matches is not None:
                inv_matches = sva_matches

            # add both matches
            if inv_matches is not None:
                # forward
                vcf_ins['me_match'] = inv_matches[0]
                ins_d = add_insertion(vcf_ins, fasta_files, alu_fasta_fh, line_fasta_fh)
                MEIs.append(ins_d)
                if 'fasta_id' in ins_d:
                    MEIs_d[ins_d['fasta_id']] = ins_d

                # reverse
                vcf_ins['me_match'] = inv_matches[1]
                ins_d = add_insertion(vcf_ins, fasta_files, alu_fasta_fh, line_fasta_fh)
                MEIs.append(ins_d)
                if 'fasta_id' in ins_d:
                    MEIs_d[ins_d['fasta_id']] = ins_d

        vcf_ins['me_match'] = me_match

        if me_match is not None:
            n_me_match += 1
            # don't expect this to happen for PAV-called L1-mediated insertions:
            if tsds['before']['len'] > tsds['after']['len']:
                warn(ins_position_str(vcf_ins) + " longer TSD sequence found before insertion point")

        # add insertions with ME match (but maybe no TSD)
        if me_match is not None:
            ins_d = add_insertion(vcf_ins, fasta_files, alu_fasta_fh, line_fasta_fh)
            MEIs.append(ins_d)
            if 'fasta_id' in ins_d:
                MEIs_d[ins_d['fasta_id']] = ins_d
                
    # close output files
    for fh in (alu_fasta_fh, line_fasta_fh):
        if fh is not None:
            fh.close()

    # run CALU, LINEU 
    run_melt_calu_lineu(MEIs_d, args.melt_jar, 'CALU', args.alu_fasta)
    run_melt_calu_lineu(MEIs_d, args.melt_jar, 'LINEU', args.line_fasta)
    for mei in MEIs:
        if 'ME_family' not in mei:
            mei['ME_family'] = mei['ME']
            mei['ME_subfamily'] = mei['ME']

    # check for overlapping annotation(s)
    # TODO - check coordinates; UCSC uses 0-start, half-open. add one to start to get 1-start, closed coords
    for mei in MEIs:
        mei['overlapping_annots'] = ""
        annots = get_overlapping_annotation(annots_by_chrom, mei['chrom'], int(mei['pos']), int(mei['pos']) + 1)
        mei['overlapping_annots'] = "|".join([":".join([a[5], a[6], a[7], a[9], a[10], a[11], a[12], a[13]])  for a in annots])
        mei['sample'] = args.sample
        
    # CSV output
    csv_fh = None
    if args.csv_output is not None:
        csv_fh = open(args.csv_output, "w")
        if csv_fh is None:
            fatal("failed to write CSV output to " + args.csv_output)
        csv_fh.write(",".join(CSV_HEADERS) + "\n")
    
        # fields for CSV output
        csv_fields = ['sample', 'chrom', 'pos', 'strand', 'ME', '%ME', '%id', '%cov',
                      'insertion_seq', 'left_flank_seq', 'right_flank_seq', 'TSD_seq', 'polyX_coords', 'ME_coords', 'insertion_coords', 'alignment',
                      # MELT CALU/LINEU
                      'ME_family', 'ME_subfamily', 'ME_start', 'ME_stop', 'ME_num_diag_matches', 'ME_num_diffs', 'ME_diffs',
                      # annotations
                      'overlapping_annots',
                      # VCF genotype and haplotype info
                      'genotype',
                      'hap1_region', 'hap2_region'
                      ]
    
        for mei in MEIs:
            csv_fh.write(",".join([mei[f] if f in mei else '' for f in csv_fields]) + "\n")
        csv_fh.close()

    # summary
    info("read " + str(len(vcf_insertions)) + " insertions from " + args.vcf)
    info("found " + str(n_tsd_before) + " insertion(s) of length >= " + str(args.min_seqlen) + " with TSDs _before_ the insertion point")
    info("found " + str(n_tsd_after) + " insertion(s) of length >= " + str(args.min_seqlen) + " with TSDs _after_ the insertion point")
    info("found " + str(n_me_match) + " insertion(s) of length >= " + str(args.min_seqlen) + " with ME matches")
    
if __name__ == '__main__':
    main()

