<script setup>
import { computed, onMounted, onUnmounted, ref } from 'vue'
import ResultExplorer from './ResultExplorer.vue'

const repository = 'https://github.com/Tomorin259/benchmark_repository'
const data = ref(null), error = ref(''), search = ref(''), directoryTab = ref('datasets')
const page = ref(location.hash.slice(1).split('?')[0] || 'home')
const initialLevel = ref(new URLSearchParams(location.hash.split('?')[1]).get('level') || 'cell')
const navigation = [['home','Overview'],['results','Results'],['directory','Datasets & Methods'],['downloads','Downloads'],['documentation','Documentation']]
const filteredDatasets = computed(() => data.value.datasets.filter(d => `${d.paper_name} ${d.primary_reference}`.toLowerCase().includes(search.value.toLowerCase())))
const filteredMethods = computed(() => data.value.methods.filter(m => `${m.manuscript_name} ${m.model_or_method}`.toLowerCase().includes(search.value.toLowerCase())))
const sourceLink = path => `${repository}/blob/${data.value.source_commit}/${path}`
const size = bytes => bytes >= 1000000 ? `${(bytes/1000000).toFixed(1)} MB` : `${(bytes/1000).toFixed(1)} KB`
const methodInput = m => ['C09','C13'].includes(m.config_id) ? 'Raw unspliced and spliced counts' : m.config_id==='C18' ? 'Existing scVelo velocity estimates' : ['C19','C20'].includes(m.config_id) ? 'Upstream method velocity estimates' : 'Unspliced and spliced RNA'
const tasks = (id, field) => [['cell','Cell'],['gene','Gene'],['kinetics','Kinetics']].filter(([key]) => data.value[key].some(r => r[field] === id)).map(([,label]) => label).join(' · ')
const downloads = computed(() => [
  ['Cell-level results', data.value.downloads.tables.filter(s => s.path.includes('/cell/'))],
  ['Gene-level results', data.value.downloads.tables.filter(s => s.path.includes('/gene/'))],
  ['Kinetic-parameter recovery', data.value.downloads.tables.filter(s => s.path.includes('/kinetics_recovery/'))],
  ['Summary figures', data.value.downloads.figures],
  ['Release metadata', data.value.downloads.metadata],
])
function focusMain(){ document.getElementById('main').focus() }
function navigate() {
  let id = location.hash.slice(1).split('?')[0] || 'home'
  if (['datasets','methods'].includes(id)) { directoryTab.value = id; id = 'directory' }
  page.value = id === 'reproduce' ? 'documentation' : id
  initialLevel.value = new URLSearchParams(location.hash.split('?')[1]).get('level') || 'cell'
  search.value = ''
  window.scrollTo(0,0)
}
async function load() {
  error.value = ''
  try {
    const response = await fetch('./data/benchmark.json')
    if (!response.ok) throw new Error(`HTTP ${response.status}`)
    data.value = await response.json()
    if (data.value.schema_version !== 1) throw new Error('Unsupported data version')
  } catch (e) { console.error(e); error.value = 'Results could not be loaded. Please try again.' }
}
onMounted(() => { window.addEventListener('hashchange',navigate); navigate(); load() })
onUnmounted(() => window.removeEventListener('hashchange',navigate))
</script>

<template>
  <a class="skip" href="#main" @click.prevent="focusMain">Skip to content</a>
  <header class="site-header"><div class="header-inner">
    <a href="#home" class="brand"><img :src="'./mark.svg'" width="36" height="36" alt=""><span>RNA velocity benchmark</span></a>
    <nav aria-label="Main navigation"><a v-for="[id,label] in navigation" :key="id" :href="`#${id}`" :class="{active:page===id}" :aria-current="page===id?'page':undefined">{{ label }}</a></nav>
    <a class="code-link" :href="repository" target="_blank" rel="noopener">Code ↗</a>
  </div></header>
  <main id="main" class="container" tabindex="-1">
    <div v-if="error" class="empty panel" role="alert"><p>{{ error }}</p><button class="button" @click="load">Try again</button></div>
    <div v-else-if="!data" class="empty" role="status">Loading results…</div>
    <template v-else>
      <template v-if="page==='home'">
        <section class="overview-header"><p class="eyebrow">Cellular trajectories · Gene kinetics · Kinetic parameters</p><h1>How well do RNA velocity methods recover biological dynamics?</h1>
          <p>Explore 20 method configurations across 14 real and three simulated datasets. The benchmark covers cell-level trajectory inference, gene-level kinetics, and kinetic-parameter recovery against scEU-seq reference values.</p>
          <div class="overview-actions"><a class="button" href="#results">Explore results</a><a class="button secondary" href="#downloads">Download results</a></div>
        </section>
        <section class="findings"><h2>Benchmark scope</h2><div class="cards">
          <article><h3>Cell-level trajectory inference</h3><p>Fourteen real datasets evaluated with five cell-level metrics.</p></article>
          <article><h3>Gene-level kinetics</h3><p>Selected real and simulated datasets evaluated with u–s phase portraits and gene-level scores.</p></article>
          <article><h3>Kinetic-parameter recovery</h3><p>Two scEU-seq metabolic-labeling datasets with experimental reference values.</p></article>
        </div></section>
        <section><h2>Browse the three evaluation tasks</h2><div class="cards overview-results">
          <a class="question-card" href="#results?level=cell"><h3>Cell-level trajectories</h3><p>Compare five cell-level metrics and inspect method-specific velocity streamlines.</p><span class="text-link">14 real datasets · 20 configurations →</span></a>
          <a class="question-card" href="#results?level=gene"><h3>Gene-level kinetics</h3><p>Inspect u–s phase portraits and gene-level scores for MURK and branching genes.</p><span class="text-link">Real and simulated settings →</span></a>
          <a class="question-card" href="#results?level=kinetics"><h3>Kinetic-parameter recovery</h3><p>Compare inferred kinetic parameters with scEU-seq experimental reference values.</p><span class="text-link">RPE1 and intestinal organoid →</span></a>
        </div></section>
        <section class="framework panel"><div><h2>Benchmarking framework</h2><p>Two evaluation levels, three complementary tasks. Each task retains its own input requirements and evaluation scope.</p></div><img :src="'./overview.png'" alt="Figure 1: RNA velocity benchmarking workflow and three evaluation tasks" loading="lazy"></section>
      </template>
      <ResultExplorer v-else-if="page==='results'" :data="data" :initial-level="initialLevel"/>
      <template v-else-if="page==='directory'">
        <div class="page-heading"><h1>Datasets & Methods</h1></div>
        <div class="tabs"><button v-for="id in ['datasets','methods']" :key="id" :class="{active:directoryTab===id}" @click="directoryTab=id; search=''">{{ id==='datasets'?'Datasets':'Methods' }}</button></div>
        <input class="directory-search" v-model="search" :aria-label="`Search ${directoryTab}`" :placeholder="`Search ${directoryTab}…`">
        <div v-if="directoryTab==='datasets'" class="directory-grid"><article v-for="d in filteredDatasets" :key="d.dataset_id" class="panel directory-card"><h2>{{ d.paper_name }}</h2><p>{{ data.dataset_descriptions[d.dataset_id] }}</p><p class="small"><strong>Available tasks:</strong> {{ tasks(d.dataset_id,'dataset') }}</p><details><summary>Dataset details</summary><p>{{ Number(d.n_obs).toLocaleString() }} cells · {{ Number(d.n_vars).toLocaleString() }} genes in the registered benchmark input.</p><p>{{ d.primary_reference }}</p><a v-if="d.authoritative_source_url.startsWith('https://')" :href="d.authoritative_source_url" target="_blank" rel="noopener">Original data source ↗</a></details><a :href="`#results?level=cell&dataset=${d.dataset_id}`" v-if="data.cell.some(r=>r.dataset===d.dataset_id)">Browse results →</a><a v-else :href="`#results?level=gene&dataset=${d.dataset_id}`">Browse results →</a></article></div>
        <div v-else class="method-grid"><article v-for="m in filteredMethods" :key="m.config_id" class="panel method-card"><h2>{{ m.manuscript_name }}</h2><p>{{ m.model_or_method }}</p><p class="small"><strong>Benchmark input:</strong> {{ methodInput(m) }}</p><p class="small"><strong>Evaluated tasks:</strong> {{ tasks(m.config_id,'method') }}</p><details><summary>Implementation details</summary><dl><dt>Version</dt><dd>{{ m.version }}</dd><dt>Execution device</dt><dd>{{ m.hardware }}</dd><dt>Configuration</dt><dd>{{ m.options }}</dd></dl></details></article></div>
        <p class="reference-link"><a :href="sourceLink('methods/README.md')" target="_blank" rel="noopener">Method execution guide ↗</a> · <a :href="sourceLink('ensemble/README.md')" target="_blank" rel="noopener">Post-hoc ensemble comparison ↗</a></p>
      </template>
      <template v-else-if="page==='downloads'">
        <div class="page-heading"><h1>Downloads</h1><p class="muted">Published score tables and summary figures. Single-method and single-gene figures can be explored in Results.</p></div>
        <section v-for="[label,files] in downloads" :key="label" class="panel downloads-section"><h2>{{ label }}</h2><div v-for="s in files" :key="s.path" class="download-row"><div><strong>{{ s.label }}</strong><p>{{ s.description }}</p><span class="small muted">{{ s.path.split('.').at(-1).toUpperCase() }} · {{ size(s.bytes) }}</span></div><a class="button secondary" :href="`./${s.download}`" download>Download</a></div></section>
        <section class="panel downloads-section"><h2>Input data & full method outputs</h2><p v-if="!data.h5ad.some(a=>a.uri)">Public H5AD downloads are being prepared. The score tables and visual previews above are available now.</p><a v-for="a in data.h5ad.filter(a=>a.uri)" :key="a.artifact_id" :href="a.uri" target="_blank" rel="noopener">{{ a.artifact_id }} ↗</a></section>
      </template>
      <template v-else-if="page==='documentation'">
        <div class="page-heading"><h1>Documentation</h1><p class="muted">Choose an operation, then follow the owning module guide for its inputs and environment.</p></div>
        <div class="reproduce-grid">
          <article class="panel"><h2>Check published results</h2><p>Validate the Cell score tables, fixed-neighbor assets, and retained figure mappings. This uses the repository bundle and does not train methods.</p><pre>python evaluation/cell/validate_results.py</pre><a :href="sourceLink('evaluation/cell/README.md')" target="_blank" rel="noopener">Cell evaluation guide ↗</a></article>
          <article class="panel"><h2>Recompute evaluation scores</h2><p>Use frozen method outputs, reference data, and the required environments. Full H5AD inputs remain external; this is separate from validating the published tables.</p><details><summary>Kinetic-parameter scoring example</summary><pre>python evaluation/kinetics_recovery/metrics.py \
  --artifact-root /path/to/kinetics-artifacts --check-only</pre><p>The evaluator computes scores from H5AD and reference inputs without saving new outputs.</p></details><div class="guide-links"><a v-for="[path,label] in [['cell','Cell'],['gene','Gene'],['kinetics_recovery','Kinetics']]" :key="path" :href="sourceLink(`evaluation/${path}/README.md`)" target="_blank" rel="noopener">{{ label }} guide ↗</a></div></article>
          <article class="panel"><h2>Regenerate summary figures</h2><p>Render summaries from the committed tables. Gene scores and kinetic trajectories are not recalculated by these commands.</p><details><summary>Gene and extended-reference examples</summary><pre>python evaluation/gene/plot_gene_assets.py --summaries-only \
  --output-root outputs/gene-summary

python -m evaluation.kinetics_recovery.expanded_reference \
  --output outputs/expanded-kinetics-summary</pre></details><a :href="sourceLink('evaluation/gene/README.md')" target="_blank" rel="noopener">Plotting guide ↗</a></article>
          <article class="panel"><h2>Evaluation scope and interpretation</h2><p>Figure 6 uses selected example genes. Figure 7 uses the canonical finite-gene summaries. The extended reference-gene analysis is a separate comparison; its candidate sets and method coverage differ.</p><p>Missing or undefined scores are excluded from means. Organoid branch means are weighted equally. The scEU-seq Dynamo comparator uses labeling information and is distinct from the RNA-seq configuration.</p><a :href="sourceLink('evaluation/kinetics_recovery/README.md')" target="_blank" rel="noopener">Read the full evaluation definitions ↗</a></article>
          <article class="panel" id="sources"><h2>Reproducibility resources</h2><p>Exact source versions, file checksums, and build metadata are available for detailed verification.</p><div class="guide-links"><a href="./data/manifest.json" target="_blank" rel="noopener">Source and checksum manifest ↗</a><a :href="sourceLink('docs/environments.md')" target="_blank" rel="noopener">Environments ↗</a><a :href="sourceLink('docs/figure_map.md')" target="_blank" rel="noopener">Paper figure sources ↗</a></div></article>
        </div>
      </template>
      <div v-else class="empty"><h1>Page not found</h1><a href="#home">Return to overview</a></div>
    </template>
  </main>
  <footer><div class="footer-inner"><span>RNA velocity benchmark</span><a href="#documentation">Documentation & reproducibility</a></div></footer>
</template>
