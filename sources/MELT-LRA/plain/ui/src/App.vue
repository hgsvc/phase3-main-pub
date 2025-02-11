<script setup>
  import { ref, computed, reactive, watch } from 'vue';
  import MEI_Viewer from './components/MEI_Viewer.vue'
  import CALU_Report from './components/CALU_Report.vue'
  import { getData } from './inline_data.js'

  // retrieve inline data
  const data_files = getData();
  
  const csv_headers = ['samples', 'chrom', 'pos', 'strand', 'ME', '%ME', '%id', '%cov', 'insertion_seq', 'left_flank_seq', 'right_flank_seq', 'TSD_seq', 'polyX_coords', 'ME_coords', 'insertion_coords', 'match_string',
            'ME_family', 'ME_subfamily', 'ME_start', 'ME_stop', 'ME_num_diag_matches', 'ME_num_diffs', 'ME_diffs', 'overlapping_annots', 'genotype', 'hap1_region', 'hap2_region']

  function parse_mei(mei_str) {
    let f = mei_str.split(',')
    let nh = csv_headers.length
    let mei = {};
    for (let i = 0;i < nh; ++i) {
        if (f[i].endsWith('%')) {
            f[i] = f[i].slice(0, -1) * 1.0
        }
        mei[csv_headers[i]] = f[i];
    }
    // overlapping annotations/repeats
    let annots = mei['overlapping_annots'].split(/\|/).filter(a => a != '')
    mei['annots'] = annots
    mei['n_overlapping_annots'] = annots.length

    // unique overlapping repeat family/families
    let overlapping_rep_fams = {}
    annots.map(a => { 
        let fam = a.split(/:/)[6]
        overlapping_rep_fams[fam]++
    })
    mei['overlapping_rep_fam'] = Object.keys(overlapping_rep_fams).sort().join(",")
    if (mei['overlapping_rep_fam'] == "") mei['overlapping_rep_fam'] = 'None'

    mei['pos'] = mei['pos'] * 1.0
    mei['insertion_size'] = mei['insertion_seq'].length
    mei['TSD_length'] = mei['TSD_seq'].length

    let pxc = mei['polyX_coords'].split('-')
    mei['polyX_length'] = pxc[1] - pxc[0] >= 0 ? pxc[1] - pxc[0] + 1 : 0
    mei['key'] = mei['chrom'] + ':' + mei['pos']
    return mei
  }

//  const mei_files = data_files.map(df => { df.file });
//  const mei_urls = mei_files.map(f => { return new URL(f, import.meta.url).href })
    
  const state = reactive({
    'data_files': data_files,
    'selected_data_file': null,
    'meis': [],
    'tab': "viewer"
  })

//  watch(() => state.selected_mei_url, (newValue) => {
//    fetch(newValue)
//    .then(res => {
//      return res.text()
//    }).then(txt => {
//      state.meis = txt.split("\n").slice(1).filter(l => !l.match(/^\s*$/)).map(ms => parse_mei(ms))
//    })
//    .catch(err => {
//      console.log("caught error " + err)
//    })
//  })

  watch(() => state.selected_data_file, (df) => {
    state.meis = df.data.split("\n").slice(1).filter(l => !l.match(/^\s*$/)).map(ms => parse_mei(ms))
  })

  state.selected_data_file = data_files[0]

</script>

<template>
  <v-app class="pa-0 ma-0 mr-3">
    <v-toolbar density="compact" color="primary" app>
      <v-app-bar-nav-icon></v-app-bar-nav-icon>
      <v-toolbar-title>
        MEI Callset: 
        <v-menu><template v-slot:activator="{ props }"><v-chip size="large" v-bind="props">{{state.selected_data_file.file}}</v-chip></template><v-list><v-list-item v-for="(i, ind) in state.data_files" @click="state.selected_data_file = i"><v-list-item-title>{{i.file}}</v-list-item-title></v-list-item></v-list></v-menu>
      </v-toolbar-title>
      <template v-slot:extension>
        <v-tabs v-model="state.tab">
          <v-tab key="viewer" value="viewer">MEI Viewer</v-tab>
          <v-tab key="report" value="report">CALU/LINEU</v-tab>
        </v-tabs>
      </template>
    </v-toolbar>
    <v-main class="pa-0 ma-0" color="white" app>
      <v-container fluid class="pa-0 ma-0">
        
        <v-card-text>
          <v-window v-model="state.tab">
            <v-window-item key="viewer" value="viewer">
               <MEI_Viewer :meis="state.meis" />
            </v-window-item>
            <v-window-item key="report" value="report">
               <CALU_Report :meis="state.meis"/>
            </v-window-item>
          </v-window>
        </v-card-text>
      </v-container>
    </v-main>
  </v-app>
</template>


<style scoped>
.v-tab {
  font-size: 1.1rem;
  font-weight: bold !important;
}
</style>
