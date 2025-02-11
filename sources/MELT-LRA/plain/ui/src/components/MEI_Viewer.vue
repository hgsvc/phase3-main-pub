<script setup>

import { ref, computed, reactive, watch } from 'vue';
import MEI from './MEI.vue'
import MiniMEI from './MiniMEI.vue'
import { format } from 'd3-format'

const props = defineProps({
    meis: Array
});

const pct_format = format(".1f")

const headers = [
{ text: 'samples', value: 'samples', sortable: true, fixed: true, width: 300 },
{ text: 'chrom', value: 'chrom', sortable: true, fixed: true, width: 80 },
{ text: 'pos', value: 'pos', sortable: true, fixed: true, width: 100 },
{ text: 'strand', value: 'strand', sortable: true, fixed: true, width: 80 },
//{ text: 'genotype', value: 'genotype', sortable: true, fixed: true },
{ text: 'ME_family', value: 'ME_family', sortable: true, fixed: true, width: 80 },
{ text: 'ME_subfamily', value: 'ME_subfamily', sortable: true, fixed: true, width: 100 },
{ text: 'ins_size', value: 'insertion_size', sortable: true, fixed: true, width: 80 },
{ text: '%ME', value: '%ME', sortable: true, fixed: true, width: 60 },
{ text: '%id', value: '%id', sortable: true, fixed: true, width: 60 },
{ text: '%coverage', value: '%cov', sortable: true, fixed: true, width: 80 },
{ text: 'TSD_bp', value: 'TSD_length', sortable: true, fixed: true, width: 80},
{ text: 'polyA/T_bp', value: 'polyX_length', sortable: true, fixed: true, width: 80},
//{ text: 'overlapping repeats', value: 'n_overlapping_annots', sortable: true, fixed: true},
{ text: 'overlapping repeat family', value: 'overlapping_rep_fam', sortable: true, fixed: true},
{ text: 'MEI', value: 'MEI', fixed: true, width: 400 }
]

// EasyDataTable pagination
const dataTable = ref()
const nextPage = () => { dataTable.value.nextPage() }
const prevPage = () => { dataTable.value.prevPage() }
const gotoPage = (p) => { if (dataTable.value != null) dataTable.value.updatePage(p) }
const maxPage = computed(() => dataTable.value?.maxPaginationNumber);
const currentPage = computed(() => dataTable.value?.currentPaginationNumber);

const sortBy = []
const sortType = []
const me_types = ['ALU', 'LINE1', 'SVA']
const me_families = ['ALU', 'AluJ', 'AluS', 'AluY', 'LINE1', 'SVA']
const genotypes = ['1|1','1|0','0|1','1|.','.|1','multiple']

const state = reactive({
    meis: [],
    selected_meis: [],
    // ME counts by type
    total_mei_counts: {},
    selected_mei_counts: {},
    selected_me_types: me_types,
    // ME counts by family
    selected_me_families: me_families,
    total_family_counts: {},
    selected_family_counts: {},
    // genotypes
    selected_genotypes: genotypes,
    total_genotype_counts: {},
    selected_genotype_counts: {},
    // overlapping repeat families
    overlapping_rep_fams: [],
    selected_overlapping_rep_fams: [],
    total_overlapping_rep_counts: {},
    selected_overlapping_rep_counts: {},
    headers: headers,
    pctid_range: [0.0, 100.0],
    ins_pctcov_range: [0.0, 100.0],
    me_pctcov_range: [0.0, 100.0],
    me_ins_length_range: [0.0, 20000.0],
    tsd_length_range: [0.0, 3000.0],
    polyx_length_range: [0.0, 2000.0],
    overlapping_repeats_range: [0.0, 10.0],
    display_mode: 'table'
})

function reset_filters() {
    state.selected_me_types = me_types
    state.selected_me_families = me_families
    state.selected_genotypes = genotypes
    state.pctid_range[0] = 0.0
    state.pctid_range[1] = 100.0
    state.ins_pctcov_range[0] = 0.0
    state.ins_pctcov_range[1] = 100.0
    state.me_pctcov_range[0] = 0.0
    state.me_pctcov_range[1] = 100.0
    state.me_ins_length_range[0] = 0.0
    state.me_ins_length_range[1] = 20000.0
    state.tsd_length_range[0] = 0.0
    state.tsd_length_range[1] = 3000.0
    state.polyx_length_range[0] = 0.0
    state.polyx_length_range[1] = 2000.0
    state.overlapping_repeats_range[0] = 0.0
    state.overlapping_repeats_range[1] = 10.0
    state.display_mode = 'table'
}

function update_selected_meis() {
    gotoPage(1)
    let f_meis = []
    
    // selected by type
    let selected_by_type = {}
    let n_selected_by_type = {}
    me_types.forEach(mt => {
        n_selected_by_type[mt] = 0
        selected_by_type[mt] = false
    })
    state.selected_me_types.forEach(mt =>{
        selected_by_type[mt] = true
    })
    
    // selected by family
    let selected_by_family = {}
    let n_selected_by_family = {}
    me_families.forEach(mf => {
        n_selected_by_family[mf] = 0
        selected_by_family[mf] = false
    })
    state.selected_me_families.forEach(mf =>{
        selected_by_family[mf] = true
    })
    
    // selected by overlapping repeat family
    let selected_by_ol_rep_fam = {}
    let n_selected_by_ol_rep_fam = {}
    
    state.overlapping_rep_fams.forEach(rf => {
        n_selected_by_ol_rep_fam[rf] = 0
        selected_by_ol_rep_fam[rf] = false
    })
    state.selected_overlapping_rep_fams.forEach(rf => {
        selected_by_ol_rep_fam[rf] = true 
    })
    
    // selected by genotype
    let selected_by_genotype = {}
    let n_selected_by_genotype = {}

    genotypes.forEach(g => {
        n_selected_by_genotype[g] = 0
        selected_by_genotype[g] = false
    })
    state.selected_genotypes.forEach(g => {
        selected_by_genotype[g] = true
    })

    state.meis.forEach(m => {
        if ((m['%id'] >= state.pctid_range[0]) && (m['%id'] <= state.pctid_range[1])
        && (m['%cov'] >= state.ins_pctcov_range[0]) && (m['%cov'] <= state.ins_pctcov_range[1])
        && (m['%ME'] >= state.me_pctcov_range[0]) && (m['%ME'] <= state.me_pctcov_range[1])
        && (m['insertion_size'] >= state.me_ins_length_range[0]) && (m['insertion_size'] <= state.me_ins_length_range[1])
        && (m['TSD_length'] >= state.tsd_length_range[0]) && (m['TSD_length'] <= state.tsd_length_range[1])
        && (m['polyX_length'] >= state.polyx_length_range[0]) && (m['polyX_length'] <= state.polyx_length_range[1])
        && (m['n_overlapping_annots'] >= state.overlapping_repeats_range[0]) && (m['n_overlapping_annots'] <= state.overlapping_repeats_range[1])
        && (selected_by_type[m['ME']])
        && (selected_by_family[m['ME_family']])
        && (selected_by_ol_rep_fam[m['overlapping_rep_fam']])
        && (selected_by_genotype[m['genotype']]))
        {
            f_meis.push(m)
            n_selected_by_type[m.ME]++
            n_selected_by_family[m.ME_family]++
            n_selected_by_ol_rep_fam[m.overlapping_rep_fam]++
            n_selected_by_genotype[m.genotype]++
        }
    })
    state.selected_meis = f_meis
    state.selected_mei_counts = n_selected_by_type
    state.selected_family_counts = n_selected_by_family
    state.selected_overlapping_rep_counts = n_selected_by_ol_rep_fam
    state.selected_genotype_counts = n_selected_by_genotype
}

function load_new_data() {
    // recompute global counts, get list of overlapping hg38 repeats
    let me_type_counts = {}
    let me_family_counts = {}
    let ol_rep_counts = {}
    let genotype_counts = {}
    
    state.meis.forEach(m => {
        if (!(m.ME in me_type_counts)) me_type_counts[m.ME] = 0
        me_type_counts[m.ME] += 1
        if (!(m.ME_family in me_family_counts)) me_family_counts[m.ME_family] = 0
        me_family_counts[m.ME_family] += 1
        if (!(m.overlapping_rep_fam in ol_rep_counts)) ol_rep_counts[m.overlapping_rep_fam] = 0
        ol_rep_counts[m.overlapping_rep_fam] += 1
        if (!(m.genotype in genotype_counts)) genotype_counts[m.genotype] = 0
        genotype_counts[m.genotype] += 1
    })
    
    state.total_mei_counts = me_type_counts
    state.total_family_counts = me_family_counts

    me_families.forEach(mf => {
        if (!(mf in me_family_counts)) {
            me_family_counts[mf] = 0
        }
    })

    state.overlapping_rep_fams = Object.keys(ol_rep_counts).sort()
    state.selected_overlapping_rep_fams = Object.keys(ol_rep_counts).sort()
    state.total_overlapping_rep_counts = ol_rep_counts
    state.total_genotype_counts = genotype_counts

    reset_filters()
    update_selected_meis()
}

watch(() => props.meis, (newValue) => { state.meis = newValue })
watch(() => state.pctid_range, (newValue) => { update_selected_meis() })
watch(() => state.ins_pctcov_range, (newValue) => { update_selected_meis() })
watch(() => state.me_pctcov_range, (newValue) => { update_selected_meis() })
watch(() => state.me_ins_length_range, (newValue) => { update_selected_meis() })
watch(() => state.tsd_length_range, (newValue) => { update_selected_meis() })
watch(() => state.polyx_length_range, (newValue) => { update_selected_meis() })
watch(() => state.overlapping_repeats_range, (newValue) => { update_selected_meis() })
watch(() => state.selected_me_types, (newValue) => { update_selected_meis() })
watch(() => state.selected_me_families, (newValue) => { update_selected_meis() })
watch(() => state.selected_overlapping_rep_fams, (newValue) => { update_selected_meis() })
watch(() => state.selected_genotypes, (newValue) => { update_selected_meis() })
watch(() => state.meis, (newValue) => { load_new_data() })

function formatRatio(n1,n2) {
    if (!n2) {
      return String(n1) + "/0";
    }
    return String(n1).padStart(String(n2).length, " ") + "/" + n2 + " (" + pct_format((n1/n2) * 100.0) + "%)";
}

function getCountRatio(m) {
    const n1 = state.selected_mei_counts[m];
    const n2 = state.total_mei_counts[m];
    return formatRatio(n1, n2)
}

function getMEColor(me) {
    if (me == 'ALU') return "#1b9e77"
    if (me == 'LINE1') return "#d95f02"
    return "#7570b3"
}

function selectAllMeiFams() {
    state.selected_me_families = me_families;
}
function deselectAllMeiFams() {
    state.selected_me_families = []
}

function selectAllOverlappingRepFams() {
    state.selected_overlapping_rep_fams = state.overlapping_rep_fams;
}
function deselectAllOverlappingRepFams() {
    state.selected_overlapping_rep_fams = []
}

function selectAllGenotypes() {
    state.selected_genotypes = genotypes
}
function deselectAllGenotypes() {
    state.selected_genotypes = []
}
state.meis = props.meis

</script>

<template>
    <v-card>
        <v-card variant="outlined" class="pa-0 ma-0">
            <v-card-title>
                <v-icon large class="pr-2">mdi-tune</v-icon><span class="font-weight-medium mr-4">Filter MEIs:</span>
                <v-chip label size="large" color="black" class="font-weight-medium ml-4">{{ formatRatio(state.selected_meis.length, state.meis.length) }}</v-chip> 
                <span class="text-h6 ml-2">total </span>
                <span v-for="me_type in ['ALU', 'LINE1', 'SVA']" class="text-h6 ml-3 mr-4">
                    <v-chip label size="large" :disabled="!state.selected_mei_counts[me_type]" :color="getMEColor(me_type)" class="font-weight-medium">{{ getCountRatio(me_type) }}</v-chip>
                    {{ me_type }}</span>
                    
                </v-card-title>
                <v-container class="pa-0 ma-0 pt-2 pl-4">
                    <v-row class="pa-0 ma-0" fluid>
                        <v-col cols="12" class="pa-0 ma-0" fluid>
                            <v-container class="pa-0 ma-0">
                                <v-row v-if="false" class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        ME type(s):
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-select v-model="state.selected_me_types" :items="['ALU', 'LINE1', 'SVA']" multiple hide-details clearable variant="outlined" density="compact" class="pa-0 ma-0 pb-2"></v-select>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        ME families:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-select v-model="state.selected_me_families" :items="me_families" multiple hide-details variant="outlined" density="compact" class="pa-0 ma-0 pb-2">
                                            <template v-slot:prepend-item>
                                                <v-list-item title="Deselect All" @click="deselectAllMeiFams"></v-list-item>
                                                <v-list-item title="Select All" @click="selectAllMeiFams"></v-list-item>
                                                <v-divider class="mt-2"></v-divider>
                                            </template>
                                            <template v-slot:selection="data">
                                                {{ data.item.value }} [{{ state.selected_family_counts[data.item.value ]}}/{{ state.total_family_counts[data.item.value] }}]
                                            </template>
                                        </v-select>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        PAV genotypes:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-select v-model="state.selected_genotypes" :items="genotypes" multiple hide-details variant="outlined" density="compact" class="pa-0 ma-0 pb-2">
                                            <template v-slot:prepend-item>
                                                <v-list-item title="Deselect All" @click="deselectAllGenotypes"></v-list-item>
                                                <v-list-item title="Select All" @click="selectAllGenotypes"></v-list-item>
                                                <v-divider class="mt-2"></v-divider>
                                            </template>
                                            <template v-slot:selection="data">
                                                <span class="pr-1">{{ data.item.value }} [{{ state.selected_genotype_counts[data.item.value ]}}/{{ state.total_genotype_counts[data.item.value] }}]</span>
                                            </template>
                                        </v-select>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        
                                    </v-col>
                                </v-row>

                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        Insertion size range:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.me_ins_length_range" :min="0" :max="20000" :step="50" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        {{ state.me_ins_length_range[0] }}bp - {{ state.me_ins_length_range[1] }}bp
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        Percent identity range:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.pctid_range" :min="0" :max="100" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        {{ state.pctid_range[0] }}% - {{ state.pctid_range[1] }}%
                                    </v-col>
                                </v-row>

                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        Insertion %coverage range:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.ins_pctcov_range" :min="0" :max="100" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        {{ state.ins_pctcov_range[0] }}% - {{ state.ins_pctcov_range[1] }}%
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        Reference ME %coverage range:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.me_pctcov_range" :min="0" :max="100" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        {{ state.me_pctcov_range[0] }}% - {{ state.me_pctcov_range[1] }}%
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        TSD length range:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.tsd_length_range" :min="0" :max="3000" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        {{ state.tsd_length_range[0] }}bp - {{ state.tsd_length_range[1] }}bp
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        polyA/T length range:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.polyx_length_range" :min="0" :max="2000" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        {{ state.polyx_length_range[0] }}bp - {{ state.polyx_length_range[1] }}bp
                                    </v-col>
                                </v-row>         
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        Number of overlapping hg38 repeats:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.overlapping_repeats_range" :min="0" :max="10" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        {{ state.overlapping_repeats_range[0] }} - {{ state.overlapping_repeats_range[1] }}
                                    </v-col>
                                </v-row>  
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        Overlapping hg38 repeat family:
                                    </v-col>
                                    <v-col cols="8" class="pa-0 ma-0">
                                        <v-select v-model="state.selected_overlapping_rep_fams" :items="state.overlapping_rep_fams" multiple hide-details variant="outlined" density="compact" class="pa-0 ma-0 pb-2">
                                            <template v-slot:prepend-item>
                                                <v-list-item title="Deselect All" @click="deselectAllOverlappingRepFams"></v-list-item>
                                                <v-list-item title="Select All" @click="selectAllOverlappingRepFams"></v-list-item>
                                                <v-divider class="mt-2"></v-divider>
                                            </template>
                                            <template v-slot:selection="data">
                                                {{ data.item.value }} [{{ state.selected_overlapping_rep_counts[data.item.value ]}}/{{ state.total_overlapping_rep_counts[data.item.value] }}]
                                            </template>
                                        </v-select>
                                    </v-col>
                                    <v-col cols="2" class="pa-0 ma-0 pl-3">
                                        
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="2" class="pa-0 ma-0">
                                        Display:
                                    </v-col>
                                    <v-col cols="10" class="pa-0 ma-0">
                                        <v-radio-group v-model="state.display_mode" inline>
                                            <v-radio label="Sortable table" value="table" density="compact"></v-radio>
                                            <v-radio label="Full-size figures [much slower]" value="figures" density="compact" class="pl-2"></v-radio>
                                        </v-radio-group>
                                    </v-col>
                                </v-row>       
                            </v-container>
                            
                        </v-col>
                    </v-row>
                </v-container>
            </v-card>
            
            <!-- sortable table view -->
            <v-card v-if="state.display_mode == 'table'" class="pa-2">
                <!-- supplemental pagination controls -->
                <div class="pa-2" style="background-color: #e0e0e0;">
                    <v-btn @click="gotoPage(1)" density="compact" :disabled="maxPage == 0">Page 1</v-btn>
                    <v-btn @click="prevPage()" density="compact" prepend-icon="mdi-arrow-left-bold" class="ml-2" :disabled="maxPage == 0">previous page</v-btn>
                    <v-btn @click="nextPage()" density="compact" append-icon="mdi-arrow-right-bold" class="ml-2" :disabled="maxPage == 0">next page</v-btn>
                    <span v-if="maxPage > 0" class="px-3 font-weight-medium">Page {{ currentPage }} / {{ maxPage }}</span>
                    <span v-else>No MEIs selected</span>
                </div>
                <EasyDataTable
                ref="dataTable"
                :headers="state.headers"
                :items="state.selected_meis"
                alternating
                border-cell
                :sort-by="sortBy"
                :sort-type="sortType"
                multi-sort
                buttons-pagination
                show-index
                class="mt-1"
                >
                <template #item-MEI="item">
                    <MiniMEI :key="item.key + '-mini'" :mei="item" />
                </template>
                
                <template #expand="item">
                    <div class="px-2">
                        <MEI :key="item.key" :mei="item" />
                        <div class="text-h6">
                            <div class="hap_div pa-1 my-1 font-weight-bold">genotype</div> {{ item.genotype }} 
                            <div class="tsd_div pa-1 my-1 ml-3 font-weight-bold">TSD</div> {{ item.TSD_seq }}

                            <div v-if="item.overlapping_annots.length > 0" class="repeat_div ml-3 pa-1 my-1 mr-2 font-weight-bold">overlapping hg38 repeats</div>
                            <span v-for="(annot, anum) in item.overlapping_annots">{{ annot }}</span>
                            <br>
                          
                            <span v-if="item.ME == 'ALU' || item.ME == 'LINE1'">
                                <div class="calu_div pa-1 my-1 font-weight-bold">{{ item.ME == 'ALU' ? 'CALU' : 'LINEU'}}</div> 
                                {{item.ME_family}} {{item.ME_subfamily}} {{item.ME_start}}-{{item.ME_stop}} diag_matches={{item.ME_num_diag_matches}} num_diffs={{ item.ME_num_diffs }}  diffs={{ item.ME_diffs }}
                            </span>


                            <div v-if="(item.hap1_region != '') || (item.hap2_region != '')" class="hap1_div pa-1 ml-3 my-1 font-weight-bold">haplotype region 1</div> {{ item.hap1_region }} 
                            <div v-if="(item.hap1_region != '') || (item.hap2_region != '')" class="hap2_div pa-1 ml-3 my-1 font-weight-bold">haplotype region 2</div> {{ item.hap2_region }}<br>

                        </div>
                    </div>
                </template>           
            </EasyDataTable>
            <!-- supplemental pagination controls -->
            <div class="pa-2" style="background-color: #e0e0e0;">
                <v-btn @click="gotoPage(1)" density="compact" :disabled="maxPage == 0">Page 1</v-btn>
                <v-btn @click="prevPage()" density="compact" prepend-icon="mdi-arrow-left-bold" class="ml-2" :disabled="maxPage == 0">previous page</v-btn>
                <v-btn @click="nextPage()" density="compact" append-icon="mdi-arrow-right-bold" class="ml-2" :disabled="maxPage == 0">next page</v-btn>
                <span v-if="maxPage > 0" class="px-3 font-weight-medium">Page {{ currentPage }} / {{ maxPage }}</span>
                <span v-else>No MEIs selected</span>
            </div>
        </v-card>
        <!-- list of figures view -->
        <v-card v-else class="pa-2" style="background-color: black;">
            <div v-for="(mei, m) in state.selected_meis">
                <MEI :mei="mei" :key="mei.key" :label="(m + 1) + '/' + state.selected_meis.length" />
            </div>
        </v-card>
    </v-card>
</template>

<style scoped>
div.tsd_div {
    display: inline-block;
    background-color: #a0ffa0;
    border: 1px solid black;
}
div.calu_div {
    display: inline-block;
    background-color: #d0d0d0;
    border: 1px solid black;
}
div.repeat_div {
    display: inline-block;
    background-color: #ffd0d0;
    border: 1px solid black;
}
div.hap_div {
    display: inline-block;
    background-color: #9090ff;
    border: 1px solid black;
}
div.hap1_div {
    display: inline-block;
    background-color: #d0d0ff;
    border: 1px solid black;
}
div.hap2_div {
    display: inline-block;
    background-color: #d0ffff;
    border: 1px solid black;
}
</style>
