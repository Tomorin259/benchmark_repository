<script setup>
import { computed, nextTick, onMounted, onUnmounted, ref } from 'vue'
import ResultExplorer from './ResultExplorer.vue'
import EnsembleResults from './EnsembleResults.vue'
import Diagnostics from './Diagnostics.vue'
import ResultDownloads from './ResultDownloads.vue'

const repository = 'https://github.com/Tomorin259/benchmark_repository'
const data = ref(null), error = ref(''), search = ref(''), directoryTab = ref('datasets')
const page = ref('home'), routeKey = ref(location.hash), initialLevel = ref('cell')
const navigation = [['home','Overview'],['results','Results'],['ensemble','Ensemble'],['diagnostics','Diagnostics'],['directory','Datasets & Methods'],['downloads','Downloads']]
const filteredDatasets = computed(() => data.value.datasets.filter(d => `${d.paper_name} ${data.value.content.datasets[d.dataset_id].technology}`.toLowerCase().includes(search.value.toLowerCase())))
const datasetGroups = computed(() => [['Real datasets',filteredDatasets.value.filter(d=>!d.dataset_id.startsWith('sim_'))],['Simulations',filteredDatasets.value.filter(d=>d.dataset_id.startsWith('sim_'))]])
const filteredMethods = computed(() => data.value.methods.filter(m => `${m.manuscript_name} ${m.model_or_method}`.toLowerCase().includes(search.value.toLowerCase())))
const sourceLink = path => `${repository}/blob/${data.value.source_commit}/${path}`
const size = bytes => bytes >= 1000000 ? `${(bytes/1000000).toFixed(1)} MB` : `${(bytes/1000).toFixed(1)} KB`
const number = value => Number(value).toLocaleString('en-US')
const rate = name => ({alpha:'α · Transcription',beta:'β · Splicing',gamma:'γ · Degradation'}[name])
const methodInput = m => ['C09','C13'].includes(m.config_id) ? 'Raw unspliced and spliced counts' : m.config_id === 'C10' ? 'Total RNA expression; native TFvelo processing' : m.config_id==='C18' ? 'Existing scVelo velocity estimates' : ['C19','C20'].includes(m.config_id) ? 'Upstream method velocity estimates' : 'Unspliced and spliced RNA'
const taskNames = {cell:'Cell',gene:'Gene',kinetics:'Kinetics'}
const taskKeys = (id, field) => Object.keys(taskNames).filter(key => (key==='cell'?[...data.value.cell,...data.value.cell_simulation]:data.value[key]).some(r => r[field] === id))
const parameterScores = id => [...new Set(data.value.kinetics.filter(r=>r.method===id).map(r=>r.parameter))].map(rate).join(' · ')
const diagnostic = id => Object.keys(data.value.diagnostics).find(key=>data.value.diagnostics[key].method===id)
const technicalFiles = computed(()=>[...data.value.downloads.metadata,...data.value.downloads.audit])
const guides = [['evaluation/cell/README.md','Cell'],['evaluation/gene/README.md','Gene'],['evaluation/kinetics_recovery/README.md','Kinetics'],['ensemble/README.md','Ensemble'],['supplement/veloae/README.md','veloAE'],['supplement/velovi/README.md','veloVI'],['supplement/tivelo_gene/README.md','TIvelo'],['supplement/multivelo/README.md','MultiVelo']]
function focusMain(){document.getElementById('main').focus()}
async function navigate(){
  const id=location.hash.slice(1).split('?')[0]||'home'
  if(['datasets','methods','simulation'].includes(id))directoryTab.value=id==='methods'?'methods':'datasets'
  page.value=['datasets','methods','simulation'].includes(id)?'directory':id==='reproduce'?'documentation':id
  initialLevel.value=new URLSearchParams(location.hash.split('?')[1]).get('level')||'cell'
  routeKey.value=location.hash;search.value='';window.scrollTo(0,0)
  await nextTick()
  if(id==='simulation')document.getElementById('simulation')?.scrollIntoView({block:'start'})
}
async function load(){
  error.value=''
  try{
    const response=await fetch('./data/benchmark.json')
    if(!response.ok)throw Error(`HTTP ${response.status}`)
    const payload=await response.json()
    if(payload.schema_version!==1||payload.content?.schema_version!==1||!payload.diagnostics)throw Error('Unsupported website data version')
    if(!Array.isArray(payload.result_packages)||!payload.result_packages.length)throw Error('Missing native result download index')
    data.value=payload;await navigate()
  }catch(e){console.error(e);error.value='Results could not be loaded. Please try again.'}
}
onMounted(()=>{window.addEventListener('hashchange',navigate);navigate();load()})
onUnmounted(()=>window.removeEventListener('hashchange',navigate))
</script>

<template>
  <a class="skip" href="#main" @click.prevent="focusMain">Skip to content</a>
  <header class="site-header"><div class="header-inner">
    <a href="#home" class="brand"><img :src="'./mark.svg'" width="36" height="36" alt=""><span>RNA velocity benchmark</span></a>
    <nav aria-label="Main navigation"><a v-for="[id,label] in navigation" :key="id" :href="`#${id}`" :class="{active:page===id}" :aria-current="page===id?'page':undefined">{{label}}</a></nav>
    <a class="code-link" :href="repository" target="_blank" rel="noopener">Code ↗</a>
  </div></header>
  <main id="main" class="container" tabindex="-1">
    <div v-if="error" class="empty panel" role="alert"><p>{{error}}</p><button class="button" @click="load">Try again</button></div>
    <div v-else-if="!data" class="empty" role="status">Loading results…</div>
    <template v-else>
      <template v-if="page==='home'">
        <section class="overview-header"><h1>RNA velocity benchmarking</h1><p>20 method configurations · 14 real datasets · 3 simulations</p></section>
        <div class="guide-links overview-links"><a href="#results?level=overall">Benchmark summary →</a><a href="#results?level=cross">Cell / Gene comparisons →</a></div>
        <section class="cards overview-results" aria-label="Evaluation tasks"><a v-for="task in data.content.scope.tasks" :key="task.id" class="question-card" :href="`#results?level=${task.id}`"><h2>{{task.label}}</h2><p>{{task.object}}</p><span class="text-link">Explore results →</span></a></section>
        <section class="framework panel"><h2>Benchmarking framework</h2><img :src="'./overview.png'" alt="RNA velocity benchmarking framework" loading="lazy"></section>
      </template>
      <ResultExplorer v-else-if="page==='results'" :key="routeKey" :data="data" :initial-level="initialLevel"/>
      <EnsembleResults v-else-if="page==='ensemble'" :key="routeKey" :data="data"/>
      <Diagnostics v-else-if="page==='diagnostics'" :key="routeKey" :data="data"/>
      <template v-else-if="page==='directory'">
        <div class="page-heading"><h1>Datasets & Methods</h1></div>
        <div class="tabs"><button v-for="id in ['datasets','methods']" :key="id" :class="{active:directoryTab===id}" :aria-pressed="directoryTab===id" @click="directoryTab=id;search=''">{{id==='datasets'?'Datasets':'Methods'}}</button></div>
        <input class="directory-search" v-model="search" :aria-label="`Search ${directoryTab}`" :placeholder="`Search ${directoryTab}…`">
        <template v-if="directoryTab==='datasets'"><section v-for="[label,items] in datasetGroups" :key="label" :id="label==='Simulations'?'simulation':undefined" class="directory-section"><h2>{{label}}</h2>
          <div class="catalog panel"><article v-for="d in items" :key="d.dataset_id" class="catalog-row">
            <div class="catalog-main"><h3>{{d.paper_name}}</h3><p class="small muted">{{data.content.datasets[d.dataset_id].technology}} · {{data.content.datasets[d.dataset_id].trajectory_pattern}}</p></div>
            <div class="catalog-size small">{{number(d.n_obs)}} cells<br>{{number(d.n_vars)}} genes</div>
            <div class="guide-links"><a v-for="key in taskKeys(d.dataset_id,'dataset')" :key="key" :href="`#results?level=${key}&dataset=${d.dataset_id}`">{{taskNames[key]}} →</a></div>
            <details class="catalog-details"><summary>Details</summary><div class="detail-body"><p>{{data.content.datasets[d.dataset_id].context}}</p><p>{{data.content.datasets[d.dataset_id].trajectory_annotation}}</p><p v-if="data.content.datasets[d.dataset_id].input_note">{{data.content.datasets[d.dataset_id].input_note}}</p>
              <p v-if="data.content.simulations[d.dataset_id]">{{data.content.simulations[d.dataset_id].normal_genes}} normal · {{data.content.simulations[d.dataset_id].murk_genes}} MURK · {{data.content.simulations[d.dataset_id].branching_genes}} branching genes</p>
              <p>{{d.primary_reference}} · {{d.primary_accession==='internal-generated'?'Study simulation':d.primary_accession}}</p>
              <a v-if="d.authoritative_source_url.startsWith('https://')" :href="d.authoritative_source_url" target="_blank" rel="noopener">Original data source ↗</a><a v-else :href="sourceLink('evaluation/gene/simulation/README.md')" target="_blank" rel="noopener">Simulation generator ↗</a>
            </div></details>
          </article><p v-if="!items.length" class="empty">No datasets match.</p></div>
        </section></template>
        <div v-else class="catalog panel"><article v-for="m in filteredMethods" :key="m.config_id" class="catalog-row">
          <div class="catalog-main"><h3>{{m.manuscript_name}}</h3><p class="small muted">{{m.model_or_method}}</p></div><div class="small muted">{{taskKeys(m.config_id,'method').map(k=>taskNames[k]).join(' · ')}}</div>
          <a v-if="diagnostic(m.config_id)" class="small" :href="`#diagnostics?method=${diagnostic(m.config_id)}`">Diagnostics →</a>
          <details class="catalog-details"><summary>Details</summary><dl class="detail-body"><dt>Benchmark input</dt><dd>{{methodInput(m)}}</dd><dt>α assumption</dt><dd>{{data.content.methods[m.config_id].model_parameters.alpha}}</dd><dt>β / γ assumption</dt><dd>{{data.content.methods[m.config_id].model_parameters.beta_gamma}}</dd><dt>Native du/dt table entry</dt><dd>{{data.content.methods[m.config_id].model_parameters.du_dt}}</dd><dt>Gene components</dt><dd>{{data.content.methods[m.config_id].gene_components}}</dd><dt>Scored kinetic parameters</dt><dd>{{parameterScores(m.config_id)||'—'}}</dd><dt>Cell-specific parameter outputs</dt><dd>{{data.content.methods[m.config_id].kinetic.parameters.map(rate).join(' · ')||'—'}}</dd><dt>Software</dt><dd>{{m.software}} · {{m.version}}</dd><dt>Device / configuration</dt><dd>{{m.hardware}} · {{m.options}}</dd><dt>Implementation</dt><dd>{{m.official_source}}</dd><template v-if="data.content.methods[m.config_id].workflow_note"><dt>Workflow</dt><dd>{{data.content.methods[m.config_id].workflow_note}}</dd></template></dl></details>
        </article><p v-if="!filteredMethods.length" class="empty">No methods match.</p></div>
        <div v-if="directoryTab==='methods'" class="panel labelled-comparator"><strong>{{data.content.labelled_control.name}}</strong><span class="small">Labeling-data comparator · Kinetics</span><a class="small" href="#results?level=kinetics">Results →</a><details><summary>Input</summary><p class="small">{{data.content.labelled_control.input}} · {{data.content.labelled_control.parameters}}</p></details></div>
      </template>
      <template v-else-if="page==='downloads'">
        <div class="page-heading"><h1>Downloads</h1></div>
        <section class="panel downloads-section" aria-labelledby="source-data-heading">
          <h2 id="source-data-heading">Source data</h2>
          <p class="small muted">v1.0.0 · Browse by figure or table and download individual CSV, TSV or NPZ files. Each directory explains its contents. The complete ZIP is optional.</p>
          <div class="guide-links">
            <a :href="data.source_data_root" target="_blank" rel="noopener">Index, plotting scripts and materials ↗</a>
            <a href="./data/files/data/SOURCE_DATA.md" download>Source Data guide</a>
            <a href="./data/files/data/source_data_manifest.tsv" download>Sizes and checksums</a>
          </div>
          <div v-for="group in data.source_data_groups" :key="group.id" class="download-row source-data-row">
            <div><strong>{{group.label}}</strong><span class="small muted">{{group.description}}</span><span class="small muted">{{number(group.files)}} files · {{size(group.bytes)}}</span></div>
            <a class="button secondary" :href="group.url" target="_blank" rel="noopener">Browse files ↗</a>
          </div>
          <div class="download-row source-data-row">
            <div><strong>Complete Source Data</strong><span class="small muted">Optional: all tables, plotting inputs, scripts and retained images · {{size(data.source_data_archive.bytes)}}</span></div>
            <a class="button secondary" :href="data.source_data_archive.url" target="_blank" rel="noopener">Download complete ZIP</a>
          </div>
        </section>
        <section class="panel downloads-section"><h2>Benchmark input files</h2><p class="small muted">25 original-format files for 17 datasets: 22 H5AD, one loom and two CSVs. Download the files you need and use them directly.</p>
          <div v-for="file in data.input_files" :key="file.artifact_id" class="download-row"><div><strong>{{file.label}}</strong><span class="small muted">{{file.format}} · {{size(file.bytes)}} · {{file.role.replaceAll('_',' ')}}</span><span class="small">{{file.path.split('/').at(-1)}}</span></div><a class="button secondary" :href="file.url" target="_blank" rel="noopener">Download file</a></div>
          <p class="small muted">Kinetic-reference Tables S1, S3 and S4 are obtained from the original study. <a :href="sourceLink('data/UPSTREAM_REFERENCES.md')" target="_blank" rel="noopener">Acquisition and verification instructions ↗</a></p>
        </section>
        <ResultDownloads :data="data"/>
        <details class="panel downloads-section"><summary>Technical resources <span>{{technicalFiles.length}} files</span></summary><p class="small"><a :href="repository" target="_blank" rel="noopener">Repository files ↗</a></p><div v-for="s in technicalFiles" :key="s.path" class="download-row"><div><strong>{{s.label}}</strong><span class="small muted">{{s.path.split('.').at(-1).toUpperCase()}} · {{size(s.bytes)}}</span></div><a class="button secondary" :href="`./${s.download}`" download>Download</a></div></details>
      </template>
      <template v-else-if="page==='documentation'">
        <div class="page-heading"><h1>Documentation</h1></div>
        <section class="panel content-panel"><h2>Metrics</h2><details v-for="(metrics,level) in data.content.metrics" :key="level" class="metric-guide"><summary>{{taskNames[level]}}</summary><div v-for="m in metrics" :key="m.key"><h3>{{m.label}}</h3><p class="small">{{m.definition}} {{m.aggregation}}</p></div></details></section>
        <section class="panel content-panel"><h2>Reproduce results</h2><div class="guide-links"><a v-for="[path,label] in guides" :key="path" :href="sourceLink(path)" target="_blank" rel="noopener">{{label}} ↗</a><a :href="sourceLink('evaluation/gene/simulation/README.md')" target="_blank" rel="noopener">Simulation ↗</a></div></section>
        <section class="panel content-panel"><h2>Data and code</h2><div class="guide-links"><a href="#downloads">Downloads →</a><a :href="sourceLink('docs/reproducibility.md')" target="_blank" rel="noopener">Reproduction guide ↗</a><a :href="sourceLink('docs/environments.md')" target="_blank" rel="noopener">Environments ↗</a><a :href="sourceLink('docs/figure_map.md')" target="_blank" rel="noopener">Figure sources ↗</a><a :href="repository" target="_blank" rel="noopener">Repository ↗</a></div></section>
      </template>
      <div v-else class="empty"><h1>Page not found</h1><a href="#home">Return to overview</a></div>
    </template>
  </main>
  <footer><div class="footer-inner"><span>RNA velocity benchmark</span><a href="#documentation">Documentation & reproducibility</a></div></footer>
</template>

<style scoped>
.source-data-row>div{min-width:0}.source-data-row strong{overflow-wrap:anywhere}.source-data-row>a{flex-shrink:0}
@media(max-width:760px){.source-data-row{flex-wrap:wrap}.source-data-row>div{flex-basis:100%}}
</style>
