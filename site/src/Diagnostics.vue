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
const methodNames = Object.fromEntries(props.data.methods.map(m=>[m.config_id,m.manuscript_name]))
const methodName = id=>methodNames[props.data.diagnostics[id].method]
const mvSourceDirectory = 'supplementary_figures/Supplementary_Figure_15/'
const mvFiles = Object.fromEntries(props.data.source_data_files.filter(r=>r.path.startsWith(mvSourceDirectory)).map(r=>[r.path.split('/').at(-1),r.url]))
const aePanel = computed(()=>props.data.diagnostics.veloae.panels.find(r=>r.dataset===dataset.value))
const format = value=>value==null?'—':Number(value).toFixed(4)
const viRows = computed(()=>props.data.diagnostics.velovi.rows.map(row=>{
  if(metric.value==='transcosine')return {dataset:row.dataset,a:row.A_native_Ms_V,b:row.B_reference_Ms_V,c:row.C_reference_Ms_affine_V,delta:row.C_minus_B}
  const values=props.data.diagnostics.velovi.conditions.map(condition=>props.data.diagnostics.velovi.full.find(r=>r.dataset===row.dataset&&r.condition===condition).highdim_cbdir)
  return {dataset:row.dataset,a:values[0],b:values[1],c:values[2],delta:values[2]-values[1]}
}))
const viMacro = computed(()=>{const r=props.data.diagnostics.velovi.macro;return {a:r.A_native_Ms_V,b:r.B_reference_Ms_V,c:r.C_reference_Ms_affine_V,delta:r.C_minus_B}})
const titles = {veloae:'Representation concentration',velovi:'Gene-wise coordinate sensitivity',tivelo:'Root-dependent orientation',multivelo:'Phase-assignment constraints'}
const guide = computed(()=>`https://github.com/Tomorin259/benchmark_repository/blob/${props.data.source_commit}/supplement/${method.value==='tivelo'?'tivelo_gene':method.value}/README.md`)
const shareUrl = computed(()=>{
  const p=new URLSearchParams({method:method.value})
  if(['veloae','velovi'].includes(method.value))p.set('dataset',dataset.value)
  else if(method.value==='tivelo'){p.set('view',view.value);if(view.value==='bonemarrow')p.set('dataset','bonemarrow')}
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
  <template v-else-if="method==='tivelo'">
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
  <template v-else-if="method==='multivelo'">
    <p class="small muted">Gate open retains MultiVelo's gene-level induction-only classification and permits repression-phase candidates during refitting. In these branching-gene comparisons, the fitting objective improved while agreement with ground-truth velocity directions decreased.</p>
    <section class="panel content-panel">
      <div class="image-heading"><h3>Representative phase portraits · {{analysis.case.gene}}</h3><a class="small" :href="mvFiles['BoostALL_5_phase.npz']" target="_blank" rel="noopener">Download phase data</a></div>
      <p class="small muted">{{names[analysis.case.dataset]}} · Supplementary Fig. 15a · CosSim labels summarize all {{analysis.case.n_cells.toLocaleString('en-US')}} cells.</p>
      <div class="multivelo-phases">
        <figure v-for="panel in analysis.phase_panels" :key="panel.label">
          <figcaption>{{panel.label}}</figcaption>
          <button class="phase-image" :aria-label="`Enlarge ${analysis.case.gene} ${panel.label} phase portrait`" @click="show(panel,`${analysis.case.gene} · ${panel.label}`)"><img :src="`./${panel.path}`" :alt="`${analysis.case.gene} spliced–unspliced phase portrait: ${panel.label}`" loading="lazy"></button>
        </figure>
      </div>
      <img class="progression-key" :src="`./${analysis.progression_key.path}`" alt="Simulated progression from early to late" loading="lazy">
      <p class="small muted compact-note">Shared spliced–unspliced coordinates and sampled cells; arrow lengths are scaled separately for display.</p>
    </section>
    <section class="panel content-panel">
      <div class="image-heading"><h3>Paired gene-level CosSim</h3><a class="small" :href="mvFiles['B_gene_scores.csv']" target="_blank" rel="noopener">Download gene scores</a></div>
      <p class="small muted">Supplementary Fig. 15b · 20 branching genes per scenario. Lines connect the same gene; filled and open points show default and modified conditions.</p>
      <button class="phase-image multivelo-pairs" aria-label="Enlarge scVelo and MultiVelo paired gene scores" @click="show(analysis.paired_scores,'scVelo and MultiVelo · Paired gene scores')"><img :src="`./${analysis.paired_scores.path}`" alt="Paired gene-level CosSim for scVelo induction-only and MultiVelo gate-open controls in normal + branching and hybrid simulations" loading="lazy"></button>
    </section>
    <section class="panel content-panel">
      <div class="image-heading"><h3>Four-condition summary</h3><a class="small" :href="mvFiles['B_condition_summary.csv']" target="_blank" rel="noopener">Download summary CSV</a></div>
      <p class="small muted">Supplementary Table 9 · Means across branching genes within each scenario; n denotes genes.</p>
      <div class="table-scroll"><table><thead><tr><th>Scenario</th><th>Method</th><th>Condition</th><th class="numeric">n</th><th class="numeric">CBDir-gene</th><th class="numeric">CosSim</th></tr></thead><tbody><tr v-for="row in analysis.rows" :key="`${row.dataset}:${row.method}:${row.condition}`"><td>{{names[row.dataset]}}</td><td>{{methodNames[row.method]}}</td><td>{{row.condition}}</td><td class="numeric">{{row.n_genes}}</td><td class="numeric">{{format(row.cbdir)}}</td><td class="numeric">{{format(row.cossim)}}</td></tr></tbody></table></div>
      <p class="reference-link"><a :href="`${data.source_data_root}/${mvSourceDirectory}`" target="_blank" rel="noopener">Figure source data ↗</a></p>
    </section>
  </template>
  <p class="reference-link"><a href="#downloads">Source data →</a></p>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>

<style scoped>
.multivelo-phases{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:20px}
.multivelo-phases figure{margin:0;min-width:0}
.multivelo-phases figcaption{font-size:14px;font-weight:600;margin-bottom:10px}
.progression-key{display:block;width:210px;max-width:100%;margin:12px auto 0}
.multivelo-pairs img{max-height:440px;object-fit:contain}
@media(max-width:760px){.multivelo-phases{grid-template-columns:1fr}}
</style>
