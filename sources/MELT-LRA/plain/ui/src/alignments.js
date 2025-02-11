// convert MEI match string to alignment spans, expressed in 0-based interbase coordinates
export function getAlignmentSpans (mei) {
    const[ins_x1, ins_x2] = mei.insertion_coords.split("-").map(x => x * 1.0)
    const[me_x1, me_x2] = mei.ME_coords.split("-").map(x => x * 1.0)
    const msl = mei.match_string.length
    var spans = []

    // offsets from left side of the alignment - will add ins_x1, me_x1 to these
    var ins_o1 = 0
    var ins_o2 = 0
    var me_o1 = 0
    var me_o2 = 0
    var n_id_bp = 0
    
    function add_span() {
        if (((ins_o2 - ins_o1) != 0) && ((me_o2 - me_o1) != 0)) {
            // handle reverse strand matches
            let ins_c = null 
            if (mei.strand == '+') {
                ins_c = [ins_x1 + ins_o1 - 1, ins_x1 + ins_o2 - 1]
            } else {
                ins_c = [ins_x2 - ins_o1, ins_x2 - ins_o2]
            }
            const span = {
                'ins': ins_c, 
                'me': [me_x1 + me_o1 - 1, me_x1 + me_o2 - 1],
                'pct_id': (n_id_bp / (ins_o2 - ins_o1)) * 100.0
            }
            spans.push(span)
            n_id_bp = 0
        }
    }

    // match string contains only ^ (gap in insertion), v (gap in ME), |, and .
    for(let i = 0;i < msl; ++i) {
        if (mei.match_string[i] == '^') {
            add_span()
            ins_o2 += 1
            ins_o1 = ins_o2
            me_o1 = me_o2
        } else if (mei.match_string[i] == 'v') {
            add_span()
            me_o2 += 1
            ins_o1 = ins_o2
            me_o1 = me_o2
        } else {
            if (mei.match_string[i] == '|') n_id_bp += 1
            // start or extend span
            ins_o2 += 1
            me_o2 += 1
        }
    }
    add_span()
    return spans
}
