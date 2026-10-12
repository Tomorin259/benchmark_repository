<script setup>
import { computed, onMounted, onUnmounted, ref } from 'vue'
import ResultExplorer from './ResultExplorer.vue'

const data = ref(null)
const error = ref('')
const page = ref(location.hash.slice(1).split('?')[0] || 'home')
const initialLevel = ref(new URLSearchParams(location.hash.split('?')[1]).get('level') || 'cell')
const search = ref('')
const repository = 'https://github.com/Tomorin259/benchmark_repository'
const navigation = [['home','Overview'],['results','Results'],['datasets','Datasets'],['methods','Methods'],['downloads','Downloads'],['reproduce','Reproduce']]
const filteredDatasets = computed(() => data.value.datasets.filter(d => `${d.paper_name} ${d.primary_reference} ${d.dataset_id}`.toLowerCase().includes(search.value.toLowerCase())))
const filteredMethods = computed(() => data.value.methods.filter(m => `${m.manuscript_name} ${m.model_or_method} ${m.config_id}`.toLowerCase().includes(search.value.toLowerCase())))
const expandedCommand = 'python -m evaluation.kinetics_recovery.expanded_reference --output outputs/expanded-kinetics-summary'
const count = value => Number(value).toLocaleString('en-US')
const size = bytes => bytes >= 1048576 ? `${(bytes / 1048576).toFixed(1)} MiB` : `${(bytes / 1024).toFixed(1)} KiB`
const sourceLink = path => `${repository}/blob/${data.value.source_commit}/${path}`
const commands = {
  cell: 'python evaluation/cell/validate_results.py',
  gene: 'python evaluation/gene/plot_gene_assets.py --summaries-only \\\n  --output-root outputs/gene-summary',
  kinetics: 'python evaluation/kinetics_recovery/metrics.py \\\n  --artifact-root /path/to/kinetics-artifacts --check-only',
}
function navigate() {
  page.value = location.hash.slice(1).split('?')[0] || 'home'
  initialLevel.value = new URLSearchParams(location.hash.split('?')[1]).get('level') || 'cell'
  search.value = ''
  window.scrollTo(0,0)
}
function focusMain() { document.getElementById('main').focus() }
async function load() {
  error.value = ''
  try {
    const response = await fetch('./data/benchmark.json')
    if (!response.ok) throw new Error(`HTTP ${response.status}`)
    const payload = await response.json()
    if (payload.schema_version !== 1) throw new Error('Unsupported data version')
    data.value = payload
  } catch (e) { error.value = `The result bundle could not be loaded (${e.message}). Please reload or rebuild the site data.` }
}
onMounted(() => { window.addEventListener('hashchange', navigate); load() })
onUnmounted(() => window.removeEventListener('hashchange', navigate))
</script>

<template>
  <a class="skip" href="#main" @click.prevent="focusMain">Skip to content</a>
  <header class="site-header"><div class="header-inner">
    <a href="#home" class="brand"><img :src="'./mark.svg'" width="36" height="36" alt=""><span>RNA velocity benchmark</span></a>
    <nav aria-label="Main navigation"><a v-for="[id,label] in navigation" :key="id" :href="`#${id}`" :class="{ active: page === id }" :aria-current="page === id ? 'page' : undefined">{{ label }}</a></nav>
    <a class="code-link" :href="repository" target="_blank" rel="noopener">Code ↗</a>
  </div></header>
  <main id="main" class="container" tabindex="-1">
    <div v-if="error" role="alert" class="panel empty"><p>{{ error }}</p><button class="button" @click="load">Try again</button></div>
    <div v-else-if="!data" class="empty" role="status">Loading benchmark results…</div>
    <template v-else>
      <template v-if="page === 'home'">
        <section class="overview-header">
          <h1>RNA velocity benchmark</h1>
          <p>20 RNA velocity method configurations across 17 datasets.</p>
          <div class="overview-actions"><a class="button" href="#results">Results</a><a class="button secondary" href="#downloads">Downloads</a><a class="button secondary" href="#reproduce">Reproduce</a></div>
        </section>
        <section class="cards overview-results">
          <a href="#results" class="question-card"><h2>Cell-level trajectory inference</h2><div class="metric-tags"><span>Consistency</span><span>ICCoh</span><span>CBDir</span><span>TransCosine</span><span>BATC</span></div></a>
          <a href="#results?level=gene" class="question-card"><h2>Gene-level kinetics evaluation</h2><div class="metric-tags"><span>CBDir-gene</span><span>CosSim</span></div></a>
          <a href="#results?level=kinetics" class="question-card"><h2>Kinetic-parameter recovery</h2><div class="metric-tags"><span>Transcription rate</span><span>Degradation rate</span></div></a>
        </section>
      </template>
      <ResultExplorer v-else-if="page === 'results'" :data="data" :initial-level="initialLevel"/>
      <template v-else-if="page === 'datasets'">
        <div class="page-heading"><h1>Datasets</h1></div>
        <input v-model="search" class="directory-search" aria-label="Search datasets" placeholder="Search datasets, references or IDs…">
        <div class="panel table-scroll"><table><thead><tr><th>Dataset</th><th>Role</th><th class="numeric">Cells</th><th class="numeric">Genes</th><th>Reference</th></tr></thead><tbody><tr v-for="d in filteredDatasets" :key="d.dataset_id"><td><strong>{{ d.paper_name }}</strong><small class="config-id">{{ d.dataset_id }}</small></td><td>{{ d.role.replaceAll('_',' ') }}</td><td class="numeric">{{ d.n_obs === '-' || !d.n_obs ? '—' : count(d.n_obs) }}</td><td class="numeric">{{ d.n_vars === '-' || !d.n_vars ? '—' : count(d.n_vars) }}</td><td><a v-if="d.authoritative_source_url.startsWith('https://')" :href="d.authoritative_source_url" target="_blank" rel="noopener">{{ d.primary_reference }} ↗</a><span v-else>{{ d.primary_reference }}</span></td></tr></tbody></table></div>

      </template>
      <template v-else-if="page === 'methods'">
        <div class="page-heading"><h1>Methods</h1></div>
        <input v-model="search" class="directory-search" aria-label="Search methods" placeholder="Search methods, models or configuration IDs…">
        <div class="method-grid"><article v-for="m in filteredMethods" :key="m.config_id" class="panel method-card"><div class="method-card-top"><span class="pill">{{ m.config_id }}</span><span class="muted small">{{ m.hardware }}</span></div><h3>{{ m.manuscript_name }}</h3><p>{{ m.model_or_method }}</p><dl><dt>Version</dt><dd>{{ m.version }}</dd><dt>Options</dt><dd>{{ m.options }}</dd></dl></article></div>
        <p class="reference-link"><a :href="sourceLink('ensemble/README.md')" target="_blank" rel="noopener">Ensemble</a></p>
      </template>
      <template v-else-if="page === 'downloads'">
        <div class="page-heading"><h1>Downloads</h1></div>
        <section v-for="[key,label] in [['tables','Result tables'],['figures','Summary figures']]" :key="key" class="panel downloads-section"><h2>{{ label }}</h2>
          <div class="table-scroll"><table><thead><tr><th>File</th><th class="numeric">{{ key === 'tables' ? 'Rows' : 'Format' }}</th><th class="numeric">Size</th><th></th></tr></thead><tbody>
            <tr v-for="s in data.downloads[key]" :key="s.path"><td>{{ s.label }}<details class="download-provenance"><summary>Source</summary><small>{{ s.path }}</small><code>{{ s.sha256 }}</code></details></td><td class="numeric">{{ key === 'tables' ? count(s.rows) : s.path.split('.').at(-1).toUpperCase() }}</td><td class="numeric">{{ size(s.bytes) }}</td><td><a :href="`./${s.download}`" download>Download ↓</a></td></tr>
          </tbody></table></div>
        </section>
        <section class="panel downloads-section"><h2>H5AD files</h2><p v-if="!data.h5ad.some(a => a.uri)" class="muted">Links pending.</p><a v-for="a in data.h5ad.filter(a => a.uri)" :key="a.artifact_id" :href="a.uri" target="_blank" rel="noopener">{{ a.artifact_id }} ↗</a></section>
      </template>
      <template v-else-if="page === 'reproduce'">
        <div class="page-heading"><h1>Reproduce</h1></div>
        <div class="reproduce-grid"><article class="panel"><h2>Cell-level results</h2><pre>{{ commands.cell }}</pre><a :href="sourceLink('evaluation/cell/README.md')" target="_blank" rel="noopener">Module guide ↗</a></article>
          <article class="panel"><h2>Gene summaries</h2><pre>{{ commands.gene }}</pre><a :href="sourceLink('evaluation/gene/README.md')" target="_blank" rel="noopener">Module guide ↗</a></article>
          <article class="panel"><h2>Kinetic parameters</h2><pre>{{ commands.kinetics }}</pre><a :href="sourceLink('evaluation/kinetics_recovery/README.md')" target="_blank" rel="noopener">Module guide ↗</a></article>
          <article class="panel"><h2>Expanded reference genes</h2><pre>{{ expandedCommand }}</pre><a :href="sourceLink('evaluation/kinetics_recovery/results/expanded_reference/README.md')" target="_blank" rel="noopener">Module guide ↗</a></article></div>
        <p class="reference-link"><a :href="sourceLink('docs/environments.md')" target="_blank" rel="noopener">Environments</a></p>
      </template>
      <div v-else class="empty"><h1>Page not found</h1><a href="#home">Return to the overview</a></div>
    </template>
  </main>
  <footer><div class="footer-inner"><span>RNA velocity benchmark</span><span v-if="data">Source <a :href="`${repository}/tree/${data.source_commit}`" target="_blank" rel="noopener">{{ data.source_commit.slice(0,7) }}</a> · <a href="./data/manifest.json">Manifest ↗</a></span></div></footer>
</template>
