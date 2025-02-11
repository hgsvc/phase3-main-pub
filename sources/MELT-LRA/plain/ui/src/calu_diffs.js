export function getCALUDiffs(mei) {
    let diffs = []
    mei.ME_diffs.split('|').forEach(d => {
        let type = null
        let start = null
        let end = null
        let seq_from = null
        let seq_to = null
        
        // case 1: deletion
        let m = d.match(/^d(\d+)(-(\d+))?$/)
        if (m) {
            type = 'd'
            start = m[1] - 1
            end = m[3] ? m[3] : m[1]
            seq_from = null
            seq_to = null
        } 
        else {
            // case 2: insertion
            m = d.match(/^i(\d+)([actg]+)$/)
            if (m) {
                type = 'i'
                start = m[1]
                end = m[1]
                seq_from = ''
                seq_to = m[2]
            }
            else {
                // case 3: substitution (always single base?)
                m = d.match(/^([actg])(\d+)([actg])$/)
                if (m) {
                    type = 's'
                    start = m[2] - 1
                    end = m[2]
                    seq_from = m[1]
                    seq_to = m[3]
                } else if ((d != "") && (d != "No Differences")) {
                    console.log("unexpected reference diff = " + d)
                }
            }
        }
        
        if (type != null) {
            let df = {
                'type': type,
                'start': start,
                'end': end,
                'seq_from': seq_from,
                'seq_to': seq_to,
                'text': d
            }
            diffs.push(df)
        }
    })
    return diffs
}
