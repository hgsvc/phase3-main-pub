#!/usr/bin/env python3

import re

ins1 = 'ACAGGGCGAGACTCCGTCTCAAAAGAAAAAAAAAAAAAAAAAAAA'
ins2 = 'TTTTTTTTTTTATTTTTAACAGGGCGAGACTCCGTCTCAAAAGAAAAAAA'
ins3 = 'ATTTTTTTTTTATTTTTAACAGGGCGAGACTCCGTCTCAAAAGAAAAAAA'
ins4 = 'TGCGCGCGCGCGCAGGGCGAGACTCCGTCTCAAAAGAAAAAATTA'
ins5 = 'AAAAAAAAAAAAAAAAAAAA'
ins6 = 'AAAAAAAAAAAAAAAAAAAT'
ins7 = 'TTTTTTTTTTTTTTTTTTTT'
ins8 = 'ATTTTTTTTTTTTTTTTTTT'

# find polyA/polyT using sliding window
def find_polyX_sw(base, start, offset, wlen, max_mismatches):
    orig_start = start
    end = start
    base_count = 0
    nonbase_count = 0
    polyX_seq = ''
    
    def add_or_remove_base(pos, add, offset):
#        print("pos=" + str(pos) + " len=" + str(len(ins)))
        nonlocal base_count
        nonlocal nonbase_count
        nonlocal polyX_seq
        
        if ins[pos] == base:
            base_count += add
        else:
            nonbase_count += add

        if add == 1:
            if offset > 0:
                polyX_seq = polyX_seq + ins[pos]
            else:
                polyX_seq = ins[pos] + polyX_seq

    if (offset == -1) and (end >= len(ins)):
        return {'len': 0, 'x1': end+2, 'x2': end+1}
        
    # start with window of size 1
    wsize = 1
    add_or_remove_base(end, 1, offset)
        
    while ((offset == 1) or (end > 0)) and ((offset == -1) or (end < len(ins) - 1)) and (nonbase_count <= max_mismatches):
#        print("start=" + str(start) + " end=" + str(end) + " wsize=" + str(wsize) + " wlen=" + str(wlen))
        # add a base to the end of the window
        end += offset
        wsize += 1

        add_or_remove_base(end, 1, offset)
            
        # remove a base from the start of the window
        if wsize > wlen:
            add_or_remove_base(start, -1, offset)
            start += offset
            wsize -= 1

    # trim nonbase characters, but only from the end (start is given and assumed to be correct?)
    nb_re = '[^' + base + ']+'
    rep = nb_re + '$' if offset > 0 else '^' + nb_re
    new_polyX_seq = re.sub(rep, '', polyX_seq, re.IGNORECASE)
    bp_trimmed = len(polyX_seq) - len(new_polyX_seq)
#    print("trimmed=" + str(new_polyX_seq) + " bp=" + str(bp_trimmed) + " start=" + str(orig_start) + " end=" + str(end))
    
    if offset > 0:
        x1 = orig_start + 1
        x2 = end + 1 - bp_trimmed
    else:
        x1 = end + 1 + bp_trimmed
        x2 = orig_start + 1

#    print("seq=" + polyX_seq)
    return { 'len': x2 - x1 + 1, 'x1': x1, 'x2': x2 }


# find polyA/polyT using exact match
def find_polyX(base, start, offset):
    end = start
    plen = 0
        
    while (end >= 0) and (end < len(ins)) and ins[end] == base:
        end += offset
        plen += 1

    if offset > 0:
        x1 = start + 1
        x2 = end 
    else:
        x1 = end + 2
        x2 = start + 1

    return { 'len': plen, 'x1': x1, 'x2': x2 }
    
def pfeat(l, f, name):
    x1 = f['x1']
    x2 = f['x2']
    fl = x2 - x1 + 1
    str = "".ljust(x1-1) + "".center(fl, "-") + "".ljust(l - x2 + 1) + "  <- " + name
    return str

# test data
for ins in (ins1, ins2, ins3, ins4, ins5, ins6, ins7, ins8):
    print("--------------------------------------------------")
    print(ins)

    # exact match
    polyA = find_polyX('A', len(ins) - 1, -1)
    polyT = find_polyX('T', 0, 1)

    # sliding window
    polyA_sw = find_polyX_sw('A', len(ins) - 1, -1, 4, 1)
    polyT_sw = find_polyX_sw('T', 0, 1, 4, 1)

    # check sliding window results against exact match
    polyA_sw_exact = find_polyX_sw('A', len(ins) - 1, -1, 1, 0)
    polyT_sw_exact = find_polyX_sw('T', 0, 1, 1, 0)

    # polyA
    if polyA['len'] > 0:
        print(pfeat(len(ins), polyA, 'polyA'))
    if polyA_sw_exact['len'] > 0:
        print(pfeat(len(ins), polyA_sw_exact, 'pA_EX'))
    if polyA_sw['len'] > 0:
        print(pfeat(len(ins), polyA_sw, 'polyA_SW'))

    # polyT
    if polyT['len'] > 0:
        print(pfeat(len(ins), polyT, 'polyT'))
    if polyT_sw_exact['len'] > 0:
        print(pfeat(len(ins), polyT_sw_exact, 'pT_EX'))
    if polyT_sw['len'] > 0:
        print(pfeat(len(ins), polyT_sw, 'polyT_SW'))

    if (polyA_sw_exact['x1'] != polyA['x1']) or (polyA_sw_exact['x2'] != polyA['x2']):
        print("** polyA mismatch")
        print("polyA=" + str(polyA))
        print("polyA_sw_exact=" + str(polyA_sw_exact))
        
    if (polyT_sw_exact['x1'] != polyT['x1']) or (polyT_sw_exact['x2'] != polyT['x2']):
        print("** polyT mismatch")
        print("polyT=" + str(polyT))
        print("polyT_sw_exact=" + str(polyT_sw_exact))

    # edge case
    polyA_last = find_polyX('A', len(ins), -1)
    print("polyA last=" + str(polyA_last))
    polyA_sw_last = find_polyX_sw('A', len(ins), -1, 1, 0)
    print("polyA SW last=" + str(polyA_sw_last))
    
    print()
