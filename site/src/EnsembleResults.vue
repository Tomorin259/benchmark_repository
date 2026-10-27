<script setup>
import { computed, nextTick, ref, watch } from 'vue'
import ImageViewer from './ImageViewer.vue'
const props = defineProps({ data: Object })
const variants = props.data.ensemble.variants
const params = new URLSearchParams(location.hash.split('?')[1])
const pairs = computed(()=>props.data.datasets.filter(d=>props.data.ensemble.scores.some(r=>r.dataset===d.dataset_id)).map(d=>({
  dataset:d, mean:props.data.ensemble.scores.find(r=>r.dataset===d.dataset_id&&r.variant===variants[0]),
  median:props.data.ensemble.scores.find(r=>r.dataset===d.dataset_id&&r.variant===variants[1]),
})))
const dataset = ref(pairs.value.some(p=>p.dataset.dataset_id===params.get('dataset'))?params.get('dataset'):'pancreas')
const rankVariant = ref(variants.includes(params.get('rank'))?params.get('rank'):variants[1])
const enlarged = ref(null), copied = ref(false), preview = ref(null)
const datasetName = computed(()=>pairs.value.find(p=>p.dataset.dataset_id===dataset.value).dataset.paper_name)
const label = variant=>variant===variants[0]?'Mean-4':'GeoMedian-4'
const format = value=>value==null?'—':Number(value).toFixed(4)
const opened = computed(()=>variants.map(variant=>({
  ...props.data.ensemble.scores.find(r=>r.dataset===dataset.value&&r.variant===variant),
  image:props.data.ensemble.streams.find(r=>r.dataset===dataset.value&&r.variant===variant),
})))
const constituentIds = props.data.ensemble.methods.map(name=>props.data.methods.find(m=>m.manuscript_name===name).config_id)
const constituentRows = computed(()=>constituentIds.map(method=>({
  ...props.data.cell.find(r=>r.dataset===dataset.value&&r.method===method),
  label:props.data.methods.find(m=>m.config_id===method).manuscript_name,
})))
const ranks = computed(()=>pairs.value.map(({dataset})=>({dataset,...Object.fromEntries(
  props.data.ensemble.rankings.filter(r=>r.dataset===dataset.dataset_id&&r.variant===rankVariant.value).map(r=>[r.metric,r]),
)})))
const shareUrl = computed(()=>`${location.origin}${location.pathname}#ensemble?${new URLSearchParams({dataset:dataset.value,rank:rankVariant.value})}`)
watch(shareUrl,value=>history.replaceState(null,'',value),{immediate:true})
async function selectDataset(id){dataset.value=id;await nextTick();preview.value.scrollIntoView({block:'start',behavior:'smooth'})}
async function copyLink(){try{await navigator.clipboard.writeText(shareUrl.value);copied.value=true;setTimeout(()=>copied.value=false,2000)}catch{copied.value=false}}
</script>

<template>
  <div class="page-heading"><h1>Post-hoc ensemble</h1><button class="button secondary" @click="copyLink">{{copied?'Copied':'Copy link'}}</button></div>
  <p class="compact-context">{{data.content.ensembles.note}}</p>
  <section ref="preview" class="analysis-preview">
    <div class="analysis-toolbar"><label for="ensemble-dataset">Dataset<select id="ensemble-dataset" v-model="dataset"><option v-for="p in pairs" :key="p.dataset.dataset_id" :value="p.dataset.dataset_id">{{p.dataset.paper_name}}</option></select></label><a :href="`#results?level=cell&dataset=${dataset}`" class="text-link">Cell results →</a></div>
    <p v-if="dataset==='forebrain'" class="small">Display: PC1/PC2 · CBDir: 30 PCs · BATC: 2D UMAP</p>
    <div class="analysis-grid"><article v-for="row in opened" :key="row.variant" class="panel analysis-image"><div class="image-heading"><h2>{{label(row.variant)}}</h2><span class="small">CBDir {{format(row.cbdir)}} · BATC {{format(row.batc)}}</span></div>
      <button class="phase-image" :aria-label="`Enlarge ${datasetName} ${label(row.variant)} streamlines`" @click="enlarged={...row.image,title:`${datasetName} · ${label(row.variant)}`} "><img :src="`./${row.image.path}`" :alt="`${datasetName} · ${label(row.variant)} streamlines`" loading="lazy"></button>
    </article></div>
    <details class="panel content-panel constituent-comparison"><summary>Compare with constituent methods</summary><div class="table-scroll"><table><thead><tr><th>Method</th><th class="numeric">CBDir</th><th class="numeric">BATC</th></tr></thead><tbody><tr v-for="row in constituentRows" :key="row.method"><td><a :href="`#results?level=cell&dataset=${dataset}&method=${row.method}`">{{row.label}} →</a></td><td class="numeric">{{format(row.cbdir)}}</td><td class="numeric">{{format(row.batc)}}</td></tr><tr v-for="row in opened" :key="row.variant" class="macro-row"><td>{{label(row.variant)}}</td><td class="numeric">{{format(row.cbdir)}}</td><td class="numeric">{{format(row.batc)}}</td></tr></tbody></table></div></details>
  </section>
  <section class="panel content-panel"><h2>Macro-average <span class="tag">14 datasets · Equal weight</span></h2><div class="table-scroll"><table>
    <thead><tr><th>Metric</th><th class="numeric">Mean-4</th><th class="numeric">GeoMedian-4</th><th class="numeric"><abbr title="GeoMedian-4 minus Mean-4">Δ</abbr></th></tr></thead>
    <tbody><tr v-for="metric in ['cbdir','batc']" :key="metric"><td>{{metric==='cbdir'?'CBDir':'BATC'}}</td><td class="numeric">{{format(data.ensemble.macro_average[variants[0]][metric])}}</td><td class="numeric">{{format(data.ensemble.macro_average[variants[1]][metric])}}</td><td class="numeric">{{format(data.ensemble.macro_average[variants[1]][metric]-data.ensemble.macro_average[variants[0]][metric])}}</td></tr></tbody>
  </table></div></section>
  <section class="panel content-panel"><h2>Dataset comparisons</h2><div class="table-scroll"><table>
    <thead><tr><th rowspan="2">Dataset</th><th class="numeric" colspan="3">CBDir</th><th class="numeric" colspan="3">BATC</th></tr><tr><template v-for="metric in ['cbdir','batc']" :key="metric"><th class="numeric">Mean-4</th><th class="numeric">GeoMedian-4</th><th class="numeric"><abbr title="GeoMedian-4 minus Mean-4">Δ</abbr></th></template></tr></thead>
    <tbody><tr v-for="p in pairs" :key="p.dataset.dataset_id" :class="{selected:p.dataset.dataset_id===dataset}"><td><button class="method-link" :aria-pressed="p.dataset.dataset_id===dataset" @click="selectDataset(p.dataset.dataset_id)">{{p.dataset.paper_name}}</button></td><template v-for="metric in ['cbdir','batc']" :key="metric"><td class="numeric">{{format(p.mean[metric])}}</td><td class="numeric">{{format(p.median[metric])}}</td><td class="numeric">{{format(p.median[metric]-p.mean[metric])}}</td></template></tr></tbody>
  </table></div></section>
  <details class="panel content-panel"><summary>Per-dataset ranks against 20 configurations</summary><label class="inline-field">Variant<select v-model="rankVariant"><option v-for="variant in variants" :key="variant" :value="variant">{{label(variant)}}</option></select></label><div class="table-scroll"><table>
    <thead><tr><th>Dataset</th><th class="numeric">CBDir rank</th><th class="numeric">BATC rank</th></tr></thead><tbody><tr v-for="r in ranks" :key="r.dataset.dataset_id"><td>{{r.dataset.paper_name}}</td><td class="numeric">{{r.cbdir.rank}} / {{r.cbdir.out_of}}<small v-if="r.cbdir.tie_count">{{r.cbdir.tie_count}} exact ties</small></td><td class="numeric">{{r.batc.rank}} / {{r.batc.out_of}}<small v-if="r.batc.tie_count">{{r.batc.tie_count}} exact ties</small></td></tr></tbody>
  </table></div></details>
  <p class="reference-link"><a href="#downloads">Source data →</a> · <a :href="`https://github.com/Tomorin259/benchmark_repository/blob/${data.source_commit}/ensemble/README.md`" target="_blank" rel="noopener">Protocol ↗</a></p>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>
