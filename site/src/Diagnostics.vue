<script setup>
import { computed, ref, watch } from 'vue'
import ImageViewer from './ImageViewer.vue'
const props = defineProps({ data: Object })
const params = new URLSearchParams(location.hash.split('?')[1])
const methods = Object.keys(props.data.diagnostics)
const method = ref(methods.includes(params.get('method'))?params.get('method'):'veloae')
const dataset = ref(props.data.diagnostics.veloae.panels.some(r=>r.dataset===params.get('dataset'))?params.get('dataset'):'pancreas')
const view = ref(params.get('view')==='bonemarrow'||params.get('dataset')==='bonemarrow'&&method.value==='tivelo'?'bonemarrow':'simulation')
const metric = ref(params.get('metric')==='highdim_cbdir'?'highdim_cbdir':'transcosine')
const enlarged = ref(null), copied = ref(false)
const analysis = computed(()=>props.data.diagnostics[method.value])
const names = Object.fromEntries(props.data.datasets.map(d=>[d.dataset_id,d.paper_name]))
const methodName = id=>props.data.methods.find(m=>m.config_id===props.data.diagnostics[id].method).manuscript_name
const aePanel = computed(()=>props.data.diagnostics.veloae.panels.find(r=>r.dataset===dataset.value))
const format = value=>value==null?'—':Number(value).toFixed(4)
const viRows = computed(()=>props.data.diagnostics.velovi.rows.map(row=>{
  if(metric.value==='transcosine')return {dataset:row.dataset,a:row.A_native_Ms_V,b:row.B_reference_Ms_V,c:row.C_reference_Ms_affine_V,delta:row.C_minus_B}
  const values=props.data.diagnostics.velovi.conditions.map(condition=>props.data.diagnostics.velovi.full.find(r=>r.dataset===row.dataset&&r.condition===condition).highdim_cbdir)
  return {dataset:row.dataset,a:values[0],b:values[1],c:values[2],delta:values[2]-values[1]}
}))
const viMacro = computed(()=>{const r=props.data.diagnostics.velovi.macro;return {a:r.A_native_Ms_V,b:r.B_reference_Ms_V,c:r.C_reference_Ms_affine_V,delta:r.C_minus_B}})
const titles = {veloae:'Representation concentration',velovi:'Gene-wise coordinate sensitivity',tivelo:'Root-dependent orientation'}
const guide = computed(()=>`https://github.com/Tomorin259/benchmark_repository/blob/${props.data.source_commit}/supplement/${method.value==='tivelo'?'tivelo_gene':method.value}/README.md`)
const shareUrl = computed(()=>{
  const p=new URLSearchParams({method:method.value})
  if(method.value!=='tivelo')p.set('dataset',dataset.value)
  else{p.set('view',view.value);if(view.value==='bonemarrow')p.set('dataset','bonemarrow')}
  if(method.value==='velovi')p.set('metric',metric.value)
  return `${location.origin}${location.pathname}#diagnostics?${p}`
})
watch(shareUrl,value=>history.replaceState(null,'',value),{immediate:true})
async function copyLink(){try{await navigator.clipboard.writeText(shareUrl.value);copied.value=true;setTimeout(()=>copied.value=false,2000)}catch{copied.value=false}}
function show(image,title){enlarged.value={...image,title}}
</script>

<template>
  <div class="page-heading"><h1>Diagnostics</h1><button class="button secondary" @click="copyLink">{{copied?'Copied':'Copy link'}}</button></div>
  <div class="tabs"><button v-for="id in methods" :key="id" :class="{active:method===id}" :aria-pressed="method===id" @click="method=id">{{methodName(id)}}</button></div>
  <div class="analysis-heading"><h2>{{titles[method]}}</h2><span class="tag">Diagnostic comparison</span><a class="small" :href="guide" target="_blank" rel="noopener">Methods ↗</a></div>
  <template v-if="method==='veloae'">
    <section class="panel content-panel"><div class="image-heading"><h3>PC1 explained variance</h3><a class="small" :href="`./${analysis.csv}`" download>Download CSV</a></div><button class="phase-image" aria-label="Enlarge PC1 explained variance" @click="show(analysis.summary,'veloAE · PC1 explained variance')"><img :src="`./${analysis.summary.path}`" alt="Raw and standardized PC1 explained variance of state and velocity across 14 datasets" class="summary-image"></button></section>
    <section class="panel content-panel"><div class="analysis-toolbar"><label for="ae-dataset">Dataset<select id="ae-dataset" v-model="dataset"><option v-for="r in analysis.panels" :key="r.dataset" :value="r.dataset">{{names[r.dataset]}}</option></select></label><button class="text-button" @click="show(analysis.plate,'veloAE · All dataset projections')">All datasets</button></div>
      <button class="phase-image diagnostic-panel-image" :aria-label="`Enlarge ${names[dataset]} PCA projection`" @click="show(aePanel,`${names[dataset]} · veloAE PCA projection`)"><img :src="`./${aePanel.path}`" :alt="`${names[dataset]} states and velocities on state-derived PC1 and PC2 axes`" loading="lazy"></button>
      <p class="small muted compact-note">Dataset-specific axes and arrow display scales.</p>
      <div class="guide-links"><a class="small" :href="`./${aePanel.path}`" download>Download image</a><a class="small" :href="`#results?level=cell&dataset=${dataset}&method=${analysis.method}`">Cell results →</a></div>
    </section>
  </template>
  <template v-else-if="method==='velovi'">
    <section class="panel content-panel"><div class="analysis-toolbar"><label for="vi-metric">Metric<select id="vi-metric" v-model="metric"><option value="transcosine">TransCosine</option><option value="highdim_cbdir">CBDir-HD</option></select></label><a class="small" :href="`./${metric==='transcosine'?analysis.csv:analysis.full_csv}`" download>Download CSV</a></div>
      <div class="condition-key"><span><strong>A</strong> Native Ms · Stored V</span><span><strong>B</strong> Reference Ms · Unchanged V</span><span><strong>C</strong> Reference Ms · Gene-wise a<sub>g</sub>V</span></div>
      <div class="table-scroll"><table><thead><tr><th>Dataset</th><th class="numeric">A</th><th class="numeric">B</th><th class="numeric">C</th><th class="numeric">C − B</th></tr></thead><tbody><tr v-for="r in viRows" :key="r.dataset" :class="{selected:r.dataset===dataset}"><td><button class="method-link" :aria-pressed="dataset===r.dataset" @click="dataset=r.dataset">{{names[r.dataset]}}</button></td><td v-for="key in ['a','b','c','delta']" :key="key" class="numeric">{{format(r[key])}}</td></tr><tr v-if="metric==='transcosine'" class="macro-row"><td>Macro-average</td><td v-for="key in ['a','b','c','delta']" :key="key" class="numeric">{{format(viMacro[key])}}</td></tr></tbody></table></div>
      <div class="guide-links"><a class="small" :href="`./${analysis.full_csv}`" download>Full results CSV</a><a class="small" :href="`#results?level=cell&dataset=${dataset}&method=${analysis.method}`">{{names[dataset]}} · Cell results →</a></div>
    </section>
  </template>
  <template v-else>
    <div class="subject-switch diagnostic-views"><button :aria-pressed="view==='simulation'" @click="view='simulation'">Simulation</button><button :aria-pressed="view==='bonemarrow'" @click="view='bonemarrow'">Bone marrow</button></div>
    <template v-if="view==='simulation'">
      <section class="panel content-panel"><div class="image-heading"><h3>Branching & Hybrid · Root and ordering</h3><a class="small" :href="`./${analysis.simulation_order.path}`" download>Download image</a></div><button class="phase-image" aria-label="Enlarge TIvelo root and ordering controls" @click="show(analysis.simulation_order,'TIvelo · Simulation root and ordering')"><img :src="`./${analysis.simulation_order.path}`" alt="Branching and hybrid roots, known progression and DPT ordering" loading="lazy"></button></section>
      <section class="panel content-panel"><div class="image-heading"><h3>Paired gene scores</h3><a class="small" :href="`./${analysis.simulation_csv}`" download>Download CSV</a></div><button class="phase-image" aria-label="Enlarge TIvelo paired gene scores" @click="show(analysis.simulation_scores,'TIvelo · Paired simulation gene scores')"><img :src="`./${analysis.simulation_scores.path}`" alt="Paired CBDir-gene and CosSim under original and known early roots" loading="lazy"></button>
        <div class="table-scroll"><table><thead><tr><th>Scenario</th><th>Refit root</th><th class="numeric">CBDir-gene</th><th class="numeric">CosSim</th><th class="numeric">Genes</th></tr></thead><tbody><tr v-for="r in analysis.simulation" :key="r.dataset+r.condition"><td>{{names[r.dataset]}}</td><td>{{r.condition==='correct_path_root'?'Original root':'Known early root'}}</td><td class="numeric">{{format(r.cbdir)}}</td><td class="numeric">{{format(r.cossim)}}</td><td class="numeric">{{r.n_genes}}</td></tr></tbody></table></div>
      </section>
    </template>
    <section v-else class="panel content-panel"><div class="image-heading"><h3>Bone marrow · Paired root refits</h3><a class="small" :href="`./${analysis.bone_csv}`" download>Download CSV</a></div><button class="phase-image" aria-label="Enlarge TIvelo bone-marrow paired root refits" @click="show(analysis.bone_image,'TIvelo · Bone-marrow paired root refits')"><img :src="`./${analysis.bone_image.path}`" alt="Mono_1 and HSC_1 root DPT on archived t-SNE coordinates and paired Cell scores" loading="lazy"></button>
      <p class="small muted compact-note">t-SNE · Map seed {{analysis.displayed_seed}} · Scores: 3 paired optimizer seeds</p>
      <div class="table-scroll"><table><thead><tr><th>Metric</th><th class="numeric">Mono_1 root</th><th class="numeric">HSC_1 root</th><th class="numeric">Paired Δ</th></tr></thead><tbody><tr v-for="r in analysis.bonemarrow" :key="r.metric"><td>{{r.metric==='cbdir'?'CBDir':'BATC'}}</td><td class="numeric">{{format(r.original_mean)}}</td><td class="numeric">{{format(r.hsc1_mean)}}</td><td class="numeric">{{format(r.delta_mean)}}</td></tr></tbody></table></div><p class="reference-link"><a :href="`#results?level=cell&dataset=bonemarrow&method=${analysis.method}`">Primary benchmark results →</a></p>
    </section>
  </template>
  <p class="reference-link"><a href="#downloads">Source data →</a></p>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>
