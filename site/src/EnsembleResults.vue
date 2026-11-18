<script setup>
import { ref } from 'vue'
import ImageViewer from './ImageViewer.vue'

const props = defineProps({ data: Object })
const selected = props.data.ensemble
const enlarged = ref(null)
const dir = './data/files/ensemble/results/selected/'
const format = value => Number(value).toFixed(4)
const delta = value => (Number(value) > 0 ? '+' : '') + format(value)
</script>

<template>
  <div class="page-heading"><h1>Five-member ensemble</h1><a class="text-link" href="#ensemble-legacy">Historical four-member results →</a></div>
  <p class="compact-context">GeoMedian-5 combines {{selected.methods.join(', ')}}. Signed neighbor weights were normalized row-wise, aggregated by the geometric median and projected with GraphVelo’s weighted-displacement formula. Member selection and evaluation used the same 14 datasets.</p>
  <section class="panel content-panel">
    <h2>Selection and performance · Figure S10</h2>
    <p>Across 14 equally weighted datasets, GeoMedian-5 scored {{format(selected.summary.cbdir)}} on CBDir and {{format(selected.summary.batc)}} on BATC.</p>
    <button class="phase-image" aria-label="Enlarge Figure S10" @click="enlarged={...selected.figure,title:'Figure S10 · Ensemble selection and performance'}"><img :src="'./'+selected.figure.path" alt="Bootstrap frequencies and dataset-level CBDir and BATC for GeoMedian-5 and 20 method configurations"></button>
    <div class="guide-links"><a :href="dir+'s10/panel_a_data.csv'" download>Screening CSV</a><a :href="dir+'s10/plot_data.csv'" download>Score CSV</a><a :href="dir+'s10/macro_means.csv'" download>Macro-average CSV</a><a v-if="data.downloads.figures.some(f=>f.path==='figures/suppFigTen.pdf')" href="./data/files/figures/suppFigTen.pdf" download>Figure S10 PDF</a></div>
    <details><summary>Candidate screening across 10,000 paired dataset resamples</summary>
      <p>Each draw identified the Pareto front of mean CBDir and BATC. The six leading candidates were separated from the remaining configurations by their front-inclusion frequencies. TSvelo was omitted from the five-member workflow on computational-cost grounds.</p>
      <div class="table-scroll"><table><thead><tr><th>Configuration</th><th class="numeric">Front inclusion</th><th>Full-data front</th><th>Member</th></tr></thead><tbody>
        <tr v-for="r in selected.screening" :key="r.config_id"><td>{{r.manuscript_name}}</td><td class="numeric">{{Number(r.bootstrap_front_percent).toFixed(2)}}%</td><td>{{r.original_mean_front==='True'?'Yes':'—'}}</td><td>{{selected.methods.includes(r.manuscript_name)?'Yes':'—'}}</td></tr>
      </tbody></table></div>
    </details>
  </section>
  <section class="panel content-panel"><h2>Matched ablations · Table S5A</h2><p class="small">Each difference is the configuration’s macro-average minus the full five-member score, calculated before rounding.</p>
    <div class="table-scroll"><table><thead><tr><th>Configuration</th><th class="numeric">CBDir</th><th class="numeric">ΔCBDir</th><th class="numeric">BATC</th><th class="numeric">ΔBATC</th></tr></thead><tbody>
      <tr v-for="r in selected.ablations" :key="r.id" :class="{'macro-row':r.id==='full5'}"><td>{{r.configuration}}</td><td class="numeric">{{format(r.cbdir)}}</td><td class="numeric">{{delta(r.delta_cbdir)}}</td><td class="numeric">{{format(r.batc)}}</td><td class="numeric">{{delta(r.delta_batc)}}</td></tr>
    </tbody></table></div><div class="guide-links"><a :href="dir+'Table_S5_A.csv'" download>Table S5A CSV</a><a :href="dir+'dataset_scores.csv'" download>All 112 dataset–configuration scores</a></div>
  </section>
  <section class="panel content-panel"><h2>Dataset scores · Table S5B</h2>
    <div class="table-scroll"><table><thead><tr><th>Dataset</th><th class="numeric">CBDir</th><th class="numeric">BATC</th></tr></thead><tbody>
      <tr v-for="r in selected.scores" :key="r.dataset"><td><a :href="'#results?level=cell&dataset='+r.dataset">{{r.display_name}}</a></td><td class="numeric">{{format(r.cbdir)}}</td><td class="numeric">{{format(r.batc)}}</td></tr>
      <tr class="macro-row"><td>Macro-average</td><td class="numeric">{{format(selected.summary.cbdir)}}</td><td class="numeric">{{format(selected.summary.batc)}}</td></tr>
    </tbody></table></div><div class="guide-links"><a :href="dir+'Table_S5_B.csv'" download>Table S5B CSV</a><a href="#downloads">Replay inputs and projected fields →</a></div>
  </section>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>
