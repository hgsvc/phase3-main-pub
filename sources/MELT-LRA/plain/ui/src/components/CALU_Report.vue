<script setup>

import { ref, computed, reactive, watch } from 'vue';
import { findOrfs } from '../orfs.js'
import { getAlignmentSpans } from '../alignments.js'
import { getCALUDiffs } from '../calu_diffs.js'
import MEI from './MEI.vue'

const props = defineProps({
    meis: Array
});

const state = reactive({
    meis: [],
    n_not_sva: null,
    n_with_deletion: null,
    n_with_overlap: null,
    n_with_subs: null,
    n_with_orphan_subs: null,
    orphan_subs_meis: []
})

// 1860 total
// 1511 Alu and LINE (254 LINE)
// 760 with deletions
// 758 with overlaps (between deletion and a S-W alignment span)
// 1481 with substitutions
// 35 with orphan subs (2.4%)

// how many have CALU/LINEU substitutions called *outside* of a S-W alignment span

function update_meis() {
    let n_not_sva = 0
    let n_with_deletion = 0
    let n_with_overlap = 0
    let n_with_subs = 0
    let n_with_orphan_subs = 0
    let orphan_subs_meis = []
    
    // search for overlaps between LINE/CALU deletions and S-W alignment spans
    state.meis.forEach(mei => {
        // me[0], me[1]
        let spans = getAlignmentSpans(mei)
        mei.n_alignments = spans.length

        // type, start, end
        let diffs = getCALUDiffs(mei)
        let deletions = diffs.filter(d => (d.type == 'd'))
        let n_deletions = deletions.length
        let has_deletions = n_deletions > 0
        
        // check for CALU/LINEU deletion overlapping with a S-W alignment span
        let n_del_overlaps = 0
        deletions.forEach(d => {
            if (d['start'] > d['end']) {
                console.log("deletion has start > end")
            }
            spans.forEach(s => {
                if (s.me[0] > s.me[1]) {
                    console.log("alignment span has start > end")
                }
                
                if ((s.me[0] < d['end']) && (s.me[1] > d['start'])) {
                    //                   console.log("alignment span = " + s.me[0] + "-" + s.me[1] + " d = " + d['start'] + " - " + d['end'])
                    n_del_overlaps += 1
                }
            })
        })
        mei.n_del_overlaps = n_del_overlaps

        // check for CALU/LINEU substitutions *not* overlapping with a S-W alignment span
        let orphan_subs = []
        let subs = diffs.filter(d => (d.type == 's'))
        let n_subs = subs.length
        let has_subs = n_subs > 0
        mei.n_subs = n_subs

        subs.forEach(sub => {
            let has_span_overlap = false
            spans.forEach(span => {
                if ((span.me[0] < sub['end']) && (span.me[1] > sub['start'])) {
                    has_span_overlap = true
                }
            })
            if (!has_span_overlap) {
                orphan_subs.push(sub)
            }
        })
        
        mei.orphan_subs = orphan_subs
        mei.n_orphan_subs = mei.orphan_subs.length

        if (n_del_overlaps > 0) {
            n_with_overlap += 1
        }
        if (mei.n_orphan_subs > 0) {
            n_with_orphan_subs += 1
            orphan_subs_meis.push(mei)
        }
        // whether the MEI has *any* CALU deletions or substitutions, respectively
        if (has_deletions) {
            n_with_deletion += 1 
        }
        if (has_subs) {
            n_with_subs += 1
        }
        if (mei.ME != 'SVA') {
            n_not_sva += 1
        }
    })
    
    state.n_not_sva = n_not_sva
    state.n_with_deletion = n_with_deletion
    state.n_with_subs = n_with_subs
    state.n_with_overlap = n_with_overlap
    state.n_with_orphan_subs = n_with_orphan_subs
    state.orphan_subs_meis = orphan_subs_meis
}

watch(() => props.meis, (newValue) => { 
    state.meis = newValue 
})
watch(() => state.meis, (newValue) => { 
    update_meis() 
})
state.meis = props.meis
        
</script>
    
    <template>
        <v-card>
            <v-card variant="outlined" class="pa-0 ma-0">
                <v-card-title>
                    {{ state.meis.length }} MEIs total<br clear="both">
                    {{ state.n_not_sva }} / {{ state.meis.length }} MEIs have CALU or LINEU calls<br clear="both">
                    {{ state.n_with_deletion }} / {{ state.n_not_sva }} MEIs have CALU/LINEU-called deletions<br clear="both">
                    {{ state.n_with_overlap }} / {{ state.n_with_deletion }} MEIs have a CALU/LINEU-called deletion that overlaps with a pipeline S-W alignment<br clear="both">
                    {{ state.n_with_subs }} / {{ state.n_not_sva }} MEIs have CALU/LINEU-called substitutions<br clear="both">
                    <br clear="both">
                    {{ state.n_with_orphan_subs }} / {{ state.n_with_subs }} MEIs have substitutions not covered by any pipeline S-W alignment:
                </v-card-title>
                <v-container class="pa-0 ma-0 pt-2 pl-4">
                    <v-row class="pa-0 ma-0">
                        <v-col cols="12" class="pa-0 ma-0">
                            
                            <div v-for="(mei, m) in state.orphan_subs_meis" class="pt-3">
                                <span style="font-weight: bold; font-size: 1rem;">{{ m+1 }} / {{  state.orphan_subs_meis.length }}</span>
                                <MEI :mei="mei" :key="mei.key" :label="(m + 1) + '/' + state.orphan_subs_meis.length" />
                                CALU/LINEU-called substitutions not covered by any pipeline S-W alignment: {{ mei.n_orphan_subs }} / {{ mei.n_subs }} {{ "[" + mei.orphan_subs.map(s => s.text).join(",") + "]" }}<br clear="both">
                                Pipeline S-W alignments that overlap with a CALU/LINEU-called deletion: {{ mei.n_del_overlaps }} / {{ mei.n_alignments }}<br clear="both">
                                {{ mei.ME_diffs }}
                            </div>
                            
                        </v-col>
                    </v-row>
                </v-container>
            </v-card>
        </v-card>
    </template> 
    
    <style scoped>
</style>
