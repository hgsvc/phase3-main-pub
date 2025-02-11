
// reverse complement NA seq
const NA_CHARS = "actgnACTGN"                                                                                                                                                                                                           
const NA_COMP = "tgacnTGACN"                                                                                                                                                                                                            
const NA_RE = new RegExp('^[' + NA_CHARS + ']*$', "i")
const NA_MAP = makemap(NA_CHARS, NA_COMP)
const START_RE = new RegExp(/^atg$/, "i")
const STOP_RE = new RegExp(/^taa|tag|tga$/, "i")

function makemap(c1, c2) {
    let cmap = {}
    let cl = c1.length
    for (let i = 0; i < cl; ++i) {
        cmap[c1[i]] = c2[i]
    }
    return cmap
}

function revcomp(seq) {
    if (!(seq.match(NA_RE))) {
        console.log("revcomp: unknown character(s) found in sequence")
    }
    let revcomp_seq =  ""
    for (let ch of seq) {
        revcomp_seq = NA_MAP[ch] + revcomp_seq
    }
    return revcomp_seq
}

function findForwardOrfs(seq) {
    let sl = seq.length
    let orfs = []
    let frames = [0,1,2]
    frames.forEach(frame =>  {
        for (let p1 = frame;p1 < sl - 2; p1 += 3) {
            let codon = seq.substring(p1, p1+3)
            let end = null;
            if (codon.match(START_RE)) {
                for (let p2 = p1 + 3;p2 < sl - 2; p2 += 3) {
                    let codon = seq.substring(p2, p2+3)
                    if (codon.match(STOP_RE)) {
                        end = p2 + 2
                        break
                    }
                }
                if (end == null) end = sl
                orfs.push({'frame': frame, 'row': frame + 1, 'start': p1 + 1, 'end': end + 1, 'len': end - p1 + 1})
                p1 = end + 1
            }
        }

    })
    return orfs
}

export function findOrfs(seq) {
    let sl = seq.length
    let fwd_orfs = findForwardOrfs(seq)
    let revseq = revcomp(seq)
    let rev_orfs = findForwardOrfs(revseq)
    // adjust coordinates of reverse-strand ORFs
    rev_orfs.forEach(o => {
      let new_start = sl - o.end + 1
      o.end = sl - o.start + 1
      o.start = new_start
      o.row += 3
    })

    return [...fwd_orfs, ...rev_orfs]
}
