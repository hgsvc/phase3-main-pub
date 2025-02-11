<script setup>

  import { reactive } from 'vue'
  import { color } from 'd3-color'
  import { interpolateYlOrRd } from 'd3-scale-chromatic'
  import { format } from 'd3-format'
  import { scaleLinear, scaleLog } from 'd3-scale'
  import { axisTop, axisBottom, axisLeft, axisRight } from 'd3-axis'
  import { findOrfs } from '../orfs.js'
  import { getAlignmentSpans } from '../alignments.js'
  import { getCALUDiffs } from '../calu_diffs.js'
  
  const props = defineProps({
    mei: { type: Object, required: true },
    label: { type: String, required: false }
  })

  const ME_LENGTHS = {
    'ALU': 281,
    'SVA': 1316,
    'SVA_A': 1387,
    'SVA_F': 1375,
    'LINE1': 6019
  }

  const mei = props.mei
  
  // SVG coordinate system
  const width = 1200
  const height = 270
  const margins = { left: 80, right: 75, top: 40, bottom: 20}
  let rx1 = margins.left
  let rx2 = width - margins.right
 
  const left_flank_bp = mei.left_flank_seq.length
  const right_flank_bp = mei.right_flank_seq.length
  const flank_bp = left_flank_bp + right_flank_bp

  // genomic coordinate system: insertion + flanking seq
  const ins_len = mei.insertion_seq.length
  const ref_xscale = scaleLinear().domain([0, ins_len + flank_bp]).range([rx1, rx2])
  
  // insertion coordinate system
  const ins_rx1 = ref_xscale(left_flank_bp)
  const ins_rx2 = ref_xscale(left_flank_bp + ins_len)
  const ins_xscale = scaleLinear().domain([0, ins_len]).range([ins_rx1, ins_rx2]);
  const ins_xaxis = axisTop(ins_xscale)
  const[px_x1, px_x2] = mei.polyX_coords.split("-").map(x => x * 1.0)
  const[ins_x1, ins_x2] = mei.insertion_coords.split("-").map(x => x * 1.0)

  // genomic flanking regions
  const lflank_xscale = scaleLinear().domain([mei.pos-left_flank_bp, mei.pos]).range([rx1, ins_rx1])
  const lflank_xaxis = axisTop(lflank_xscale).tickValues([mei.pos]).tickFormat(bp => mei.chrom + ":" + bp).tickSizeOuter(0)
  const rflank_xscale = scaleLinear().domain([mei.pos + ins_len, mei.pos + ins_len + right_flank_bp]).range([ins_rx2, rx2])
  const rflank_xaxis = axisTop(rflank_xscale).tickValues([]).tickFormat(bp => "").tickSizeOuter(0)

  // reference ME coordinate system
  const me_len = ME_LENGTHS[mei.ME]
 
  // srhink ME to fit
  let me_rx1 = ins_rx1
  let me_rx2 = ins_rx2
  
  // but use same scale if ME is smaller
  if (me_len <= ins_len + flank_bp) {
    const me_diff = (ins_len + flank_bp) - me_len
    me_rx1 = ref_xscale(me_diff/2)
    me_rx2 = ref_xscale(me_diff/2 + me_len)
  } 

  const me_domain = mei.strand == '+' ? [0, me_len] : [me_len, 0]
  const me_xscale = scaleLinear().domain(me_domain).range([me_rx1, me_rx2]);
  const me_xaxis = axisBottom(me_xscale)
  const[me_x1, me_x2] = mei.ME_coords.split("-").map(x => x * 1.0)

  // set color based on percent identity
  const color_fn = interpolateYlOrRd;
  let col = null;

  function pctid_color(pctid, alpha) {
    const clr = color(color_fn(pctid / 100.0))
    clr.opacity = alpha
    return clr.formatRgb()
  }

  const feat_y_offset = 6
  const lbl_y_offset = 22
  const ins_axis_y = margins.top + 50
  const me_axis_y = height - margins.bottom - 50

  // convert spans to polygons
  const match_y1 = ins_axis_y + feat_y_offset + 6
  const match_y2 = me_axis_y - feat_y_offset - 4

  let spans = getAlignmentSpans(mei)
  spans.forEach(s => {
    s.points = []
    s.points.push([ins_xscale(s.ins[0]), match_y1])
    s.points.push([ins_xscale(s.ins[1]), match_y1])
    s.points.push([me_xscale(s.me[1]), match_y2])
    s.points.push([me_xscale(s.me[0]), match_y2])
    s.points_str = s.points.join(" ")
  })

  // convert ME_diffs string to list of differences (with the reference)
  // e.g. "d1-137|g211a|c236t|i252gcagtcc|c248g|a262g|a252g"
  let diffs = getCALUDiffs(mei)
  
  const diffs_y1 = match_y2 + 6
  const diffs_y2 = match_y2 + 16
  const base_colors = { 'a': '#209af7', 'c': '#fa3c4c', 'g': '#58fa3c', 't': '#faed3c'}

  diffs.forEach(d => { 
    let x_start = me_xscale(d['start'])
    let x_end = me_xscale(d['end'])          
    d['y1'] = diffs_y1
    d['y2'] = diffs_y2
    d['x1'] = x_start <= x_end ? x_start : x_end
    d['x2'] = x_start <= x_end ? x_end : x_start
    if (d['type'] == 's') d['color'] = base_colors[d['seq_to']]
  })

  const me_ref_str = ": " + mei['%ME'] + "% ME coverage at " + mei['%id'] + "% identity"

  // color key - percent identity
  const cb_height = 20
  const cb_width = 15
  const cb_y = height - margins.bottom - (cb_height * 5)
  const color_key_blocks = [
    { 'id': 95, 'color': pctid_color(95, 1), 'y': cb_y, height: cb_height },
    { 'id': 90, 'color': pctid_color(90, 1), 'y': cb_y + cb_height, 'height': cb_height },
    { 'id': 80, 'color': pctid_color(80, 1), 'y': cb_y + cb_height * 2, 'height': cb_height },
    { 'id': 70, 'color': pctid_color(70, 1), 'y': cb_y + cb_height * 3, 'height': cb_height },
    { 'id': 60, 'color': pctid_color(60, 1), 'y': cb_y + cb_height * 4, 'height': cb_height },
  ]

  // color key - base substitutions
  const base_color_key_blocks = [
    { 'base': 'A', 'color': base_colors['a'], 'y': cb_y, height: cb_height },
    { 'base': 'C', 'color': base_colors['c'], 'y': cb_y + cb_height, 'height': cb_height },
    { 'base': 'G', 'color': base_colors['g'], 'y': cb_y + cb_height * 2, 'height': cb_height },
    { 'base': 'T', 'color': base_colors['t'], 'y': cb_y + cb_height * 3, 'height': cb_height },
    { 'base': 'del', 'color': '#202020', 'y': cb_y + cb_height * 4, 'height': cb_height },
  ]

  // L1s only - find ORFs
  let orfs = []
  if (mei.ME == 'LINE1') {
    orfs = findOrfs(mei['insertion_seq'])
    orfs.forEach(o => {
      o.x1 = ins_xscale(o.start)
      o.x2 = ins_xscale(o.end)
    })
  }

  const state = reactive({
    // SVG
    width: width,
    height: height,
    margins: margins,

    // props
    label: props.label,
    // MEI
    mei: mei, 
    ins_len: ins_len,
    tsd_len: mei.TSD_seq.length,
    // flanking regions
    rflank_x: mei.pos + ins_len,
    lflank_xscale: lflank_xscale,
    lflank_xaxis: lflank_xaxis,
    rflank_xscale: rflank_xscale,
    rflank_xaxis: rflank_xaxis,
    // insertion and associated features
    ins_xscale: ins_xscale, 
    ins_xaxis: ins_xaxis,
    polyx_x1: px_x1,
    polyx_x2: px_x2,
    ins_x1: ins_x1,
    ins_x2: ins_x2,
    // reference ME
    me_len: me_len,
    me_xscale: me_xscale,
    me_xaxis: me_xaxis,
    me_x1: me_x1,
    me_x2: me_x2,
    // y offsets
    // insertion feature on ref genome axis (top)
    ins_label_y: margins.top - 10,
    ins_y: margins.top + 3,
    // axes for left and right genomic flanking regions (top)
    flank_axis_y: margins.top,
    flank_feat_y: margins.top + feat_y_offset,
    flank_feat_lbl_y: margins.top + lbl_y_offset,
    // insertion axis and features (below genomic axis)
    ins_axis_y: ins_axis_y,
    ins_feat_y: ins_axis_y + feat_y_offset,
    ins_feat_lbl_y: ins_axis_y + lbl_y_offset,
    // mobile element reference 
    me_axis_y: me_axis_y,
    me_lbl_y: me_axis_y + 50,
    // alignment spans
    spans: spans,
    me_ref_str: me_ref_str,
    // differences with reference
    diffs: diffs,
    // ORFs
    orfs: orfs,
    // color keys
    color_key_blocks: color_key_blocks,
    base_color_key_blocks: base_color_key_blocks
  })
  
</script>

<template>
  <div class="greetings pa-0">
    <svg
      ref="mei_svg"
      key="mei-svg-1"
      :width="state.width"
      :height="state.height"
      class="pa-0 ma-0"
      style="border: 1px solid white; background-color: black; color: rgba(235,235,235,0.64);"
      >

      <!-- ad-hoc color key / percent identity -->
      <rect v-for="ck in color_key_blocks" :x="10" :y="ck.y" :width="20" :height="ck.height" :fill="ck.color"></rect>
      <text v-for="ck in color_key_blocks" :x="35" :y="ck.y + ck.height -5" fill="#d0d0d0">{{ck.id + "%"}}</text>

      <!-- color key for base substitutions -->
      <rect v-for="ck in base_color_key_blocks" :x="width - 30" :y="ck.y" :width="20" :height="ck.height-2" :fill="ck.color" :stroke="ck.base == 'del' ? '#ffffff' : 'none'" stroke-dasharray="2,4,2,4"></rect>
      <text v-for="ck in base_color_key_blocks" :x="width - 40" :y="ck.y + ck.height -5" fill="#d0d0d0" text-anchor="end">{{ck.base}}</text>

      <!-- ORFS -->
      <text v-if="state.orfs.length > 0" :x="10" :y="120" font-weight="bold" fill="#f189f5" stroke="#f189f5" font-size="1rem">ORFs</text>

      <!-- insertion feature on ref genome -->
      <line :x1="state.ins_xscale(0)" :x2="state.ins_xscale(state.ins_len)" :y1="state.ins_y" :y2="state.ins_y" stroke-width="3" stroke="#ffffff" />
      <text class="ins_label" :x="state.ins_xscale(state.ins_len/2)" :y="state.ins_label_y" fill="#ffffff">{{ state.ins_len + " bp insertion [" + state.mei.strand + " strand]" }}</text>
      <!-- genomic sequence flanking regions -->
      <g v-axis="state.lflank_xaxis" class="xaxis" :transform="`translate(0,${state.flank_axis_y})`"></g>
      <g v-axis="state.rflank_xaxis" class="xaxis" :transform="`translate(0,${state.flank_axis_y})`"></g>

      <!-- TSD in right flanking sequence-->
      <line v-if="state.tsd_len > 0" :x1="state.rflank_xscale(state.rflank_x)" :x2="state.rflank_xscale(state.rflank_x + state.tsd_len)" :y1="state.flank_feat_y" :y2="state.flank_feat_y" stroke-width="4" stroke="#a0ffa0"/>
      <text v-if="state.tsd_len > 0" :x="state.rflank_xscale(state.rflank_x)" :y="state.flank_feat_lbl_y" fill="#ffffff">TSD</text>

      <!-- inserted sequence axis -->
      <g v-axis="state.ins_xaxis" class="xaxis" :transform="`translate(0,${state.ins_axis_y})`">
      </g>

      <!-- reference ME axis -->
      <g v-axis="state.me_xaxis" class="xaxis" :transform="`translate(0,${state.me_axis_y})`">
      </g>
      <text v-if="state.mei.strand == '+'" class="me_label" :x="state.me_xscale(0)" :y="state.me_lbl_y" fill="#ffffff" text-anchor="start">{{ state.me_len + " bp " + state.mei.ME + " reference " + me_ref_str + " >>>"}}</text>
      <text v-else class="me_label" :x="state.me_xscale(0)" :y="state.me_lbl_y" fill="#ffffff" text-anchor="end">{{ "<<< " + state.me_len + " bp " + state.mei.ME + " reference " + me_ref_str}}</text>

      <!-- match/alignment spans -->
      <polygon v-for="s in state.spans" :points="s.points_str" :fill="pctid_color(s.pct_id, 0.5)" :stroke="pctid_color(s.pct_id, 1.0)" stroke-width="2" />

      <!-- CALU/LINEU reference diffs - deletions -->
      <rect v-for="d in state.diffs.filter(d => d.type == 'd')" :x="d.x1" :y="d.y1" :width="d.x2 - d.x1" :height="d.y2 - d.y1" fill="#202020" stroke="#ffffff" stroke-dasharray="2,4,2,4"></rect>
      <!-- CALU/LINEU reference diffs - substitutions -->
      <line v-for="d in state.diffs.filter(d => d.type == 's')" :x1="d.x1" :x2="d.x1" :y1="d.y1" :y2="d.y2" :stroke="d.color" stroke-width="2" />

      <!-- ORFS in insertion sequence -->
      <rect v-for="o in state.orfs.filter(o => (o.len >= 100))" :x="o.x1" :y="state.ins_feat_y + 2 + (o.row * 6)" :width="o.x2 - o.x1" :height="4" fill="#f189f5" stroke="#f189f5"></rect>

      <!-- TSD, polyX in ME coords -->
      <line v-if="state.tsd_len > 0" :x1="state.ins_xscale(0)" :x2="state.ins_xscale(state.tsd_len)" :y1="state.ins_feat_y" :y2="state.ins_feat_y" stroke-width="4" stroke="#8ccf9e"/>
      <text v-if="state.tsd_len > 0" :x="state.ins_xscale(0)" :y="state.ins_feat_lbl_y" fill="#ffffff">TSD</text>
      <line v-if="state.mei.polyX_length > 0" :x1="state.ins_xscale(state.polyx_x1)" :x2="state.ins_xscale(state.polyx_x2)" :y1="state.ins_feat_y" :y2="state.ins_feat_y" stroke-width="4" stroke="#8080ff"/>
      <text v-if="state.mei.polyX_length > 0" :x="state.ins_xscale(state.polyx_x1)" :y="state.ins_feat_lbl_y" fill="#ffffff">{{ state.polyx_x2 == state.ins_len ? "polyA" : "polyT" }}</text>
      </svg>
  </div>
</template>

<style scoped>
.xaxis {
  font-size: 18px;
}

.ins_label {
  font-size: 18px;
  text-anchor: middle;
}

.me_label {
  font-size: 18px;
}

</style>
