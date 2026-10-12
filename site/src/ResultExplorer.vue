<script setup>
import { computed, ref, watch } from 'vue'
import GeneVisual from './GeneVisual.vue'

const props = defineProps({ data: Object, initialLevel: { type: String, default: 'cell' } })
const mode = ref('cell')
const dataset = ref('pancreas')
const metric = ref('cbdir')
const scope = ref('expanded')
const comparison = ref('available')
const branch = ref('')
const parameter = ref('alpha')
const gene = ref('')
const query = ref('')
const order = ref('method')
const display = ref('table')
const visualIndex = ref(null)
const visualLoading = ref(false)
const visualError = ref('')
watch(display, async value => {
  if (value !== 'visual' || visualIndex.value || visualLoading.value || !props.data.visuals) return
  visualLoading.value = true
  visualError.value = ''
  try {
    const response = await fetch(`./${props.data.visuals.index}`)
    if (!response.ok) throw new Error(`HTTP ${response.status}`)
    const index = await response.json()
    if (index.schema_version !== 1) throw new Error('Unsupported visual data version')
    visualIndex.value = index
  } catch (e) { visualError.value = `Visuals could not be loaded (${e.message}).` }
  finally { visualLoading.value = false }
})
const selected = ref([...props.data.methods.map(m => m.config_id), 'dynamo_labelled'])
const names = Object.fromEntries(props.data.methods.map(m => [m.config_id, m.manuscript_name]))
names.dynamo_labelled = 'Dynamo (label data)'
const tables = computed(() => props.data[mode.value])
const datasets = computed(() => props.data.datasets.filter(d => tables.value.some(r => r.dataset === d.dataset_id)))
const datasetName = computed(() => props.data.datasets.find(d => d.dataset_id === dataset.value)?.paper_name)
const branches = computed(() => [...new Set(props.data.kinetics.filter(r => r.dataset === dataset.value && r.scope === scope.value).map(r => r.branch))].sort())
const parameters = computed(() => [...new Set(props.data.kinetics.filter(r => r.dataset === dataset.value && r.scope === scope.value && r.branch === branch.value).map(r => r.parameter))].sort())
const commonGenes = computed(() => props.data.kinetics_common_membership.find(r =>
  r.dataset === dataset.value && r.branch === branch.value && r.parameter === parameter.value && r.comparison === comparison.value)?.genes ?? [])
const commonMethods = computed(() => props.data.kinetics_summary.filter(r => r.scope === 'expanded' &&
  r.dataset === dataset.value && r.branch === branch.value && r.parameter === parameter.value && r.comparison === comparison.value).map(r => r.method))
const candidates = computed(() => tables.value.filter(r => r.dataset === dataset.value && (mode.value !== 'kinetics' ||
  (r.scope === scope.value && r.branch === branch.value && r.parameter === parameter.value &&
    (scope.value !== 'expanded' || comparison.value === 'available' ||
      (commonGenes.value.includes(r.gene) && commonMethods.value.includes(r.method)))))))
const genes = computed(() => [...new Set(candidates.value.map(r => r.gene).filter(Boolean))].sort())
const metrics = computed(() => mode.value === 'cell' ? [
  ['cbdir', 'CBDir'], ['iccoh', 'ICCoh'], ['consistency', 'Consistency'], ['transcosine', 'TransCosine'], ['batc', 'BATC'],
] : mode.value === 'gene' ? [['cbdir', 'CBDir-gene'], ...(dataset.value.startsWith('sim_') ? [['cossim', 'CosSim']] : [])] : [['score', 'Spearman correlation']])
const metricName = computed(() => metrics.value.find(m => m[0] === metric.value)?.[1])
const comparisonNames = { available: 'Available genes', core_conventional: 'Common: core conventional',
  all_conventional: 'Common: all conventional', all_displayed_methods: 'Common: all displayed methods' }

watch(mode, value => {
  dataset.value = value === 'kinetics' ? 'sceu_rpe1' : 'pancreas'
  metric.value = value === 'kinetics' ? 'score' : 'cbdir'
  query.value = ''
})
watch(() => props.initialLevel, value => {
  mode.value = ['cell', 'gene', 'kinetics'].includes(value) ? value : 'cell'
}, { immediate: true })
watch(branches, value => { if (!value.includes(branch.value)) branch.value = value[0] ?? '' }, { immediate: true })
watch(parameters, value => { if (!value.includes(parameter.value)) parameter.value = value[0] ?? 'alpha' }, { immediate: true })
watch(scope, value => { if (value !== 'expanded') comparison.value = 'available' })
watch([mode, dataset], () => {
  gene.value = mode.value !== 'cell' && scope.value !== 'selected' ? '__summary__' : ''
  display.value = 'table'
})
watch([mode, dataset, scope, branch, parameter, genes, metric], () => {
  const summaryAllowed = (mode.value === 'kinetics' && scope.value !== 'selected') || (mode.value === 'gene' && metric.value === 'cbdir')
  if (!genes.value.includes(gene.value) && !(gene.value === '__summary__' && summaryAllowed)) gene.value = summaryAllowed ? '__summary__' : genes.value[0] ?? ''
  if (!metrics.value.some(m => m[0] === metric.value)) metric.value = metrics.value[0][0]
}, { immediate: true })
watch([mode, gene], () => { if (gene.value === '__summary__' || mode.value === 'cell') display.value = 'table' })

const unfilteredRows = computed(() => {
  const data = mode.value === 'kinetics' && gene.value === '__summary__'
    ? props.data.kinetics_summary.filter(r => r.dataset === dataset.value && r.branch === branch.value &&
      r.parameter === parameter.value && r.scope === scope.value && r.comparison === comparison.value)
    : mode.value === 'gene' && gene.value === '__summary__'
      ? props.data.gene_summary.filter(r => r.dataset === dataset.value)
      : candidates.value.filter(r => mode.value === 'cell' || r.gene === gene.value)
  return data.map(r => ({ ...r, name: names[r.method], value: r[metric.value] }))
})
const availableMethods = computed(() => [...new Set(unfilteredRows.value.map(r => r.method))])
const rows = computed(() => unfilteredRows.value.filter(r => selected.value.includes(r.method) &&
  `${r.name} ${r.method}`.toLowerCase().includes(query.value.toLowerCase())).sort((a, b) => {
    if (order.value === 'method') {
      const ids = mode.value === 'kinetics' ? props.data.kinetics_method_order : props.data.methods.map(m => m.config_id)
      return ids.indexOf(a.method) - ids.indexOf(b.method)
    }
    if (a.value === null) return b.value === null ? a.name.localeCompare(b.name) : 1
    if (b.value === null) return -1
    return (order.value === 'descending' ? b.value - a.value : a.value - b.value) || a.name.localeCompare(b.name)
  }))
const finiteCount = computed(() => rows.value.filter(r => r.value !== null).length)
const scopeLabel = computed(() => mode.value === 'cell' ? 'Dataset · reference 30-NN'
  : mode.value === 'gene' ? gene.value === '__summary__' ? 'Displayed genes · finite-gene mean · reference 30-NN' : 'Gene · reference 30-NN'
  : scope.value === 'expanded' ? `Expanded reference · ${props.data.kinetics_expanded_scope[dataset.value]} candidate genes · ${comparisonNames[comparison.value]}`
  : gene.value === '__summary__' ? 'All valid genes · mean'
  : scope.value === 'selected' ? 'Selected genes · single gene' : 'All evaluated genes · single gene')
const sourcePaths = computed(() => [...new Set([...rows.value.map(r => mode.value === 'gene' && metric.value === 'cossim'
  ? 'evaluation/gene/results/simulation/gene_cbdir_cossim_reference.csv' : r.source),
  ...(mode.value === 'kinetics' && scope.value === 'expanded' ? ['evaluation/kinetics_recovery/results/expanded_reference/manifest.json',
    'evaluation/kinetics_recovery/results/expanded_reference/gene_scope.csv',
    ...(comparison.value !== 'available' ? ['evaluation/kinetics_recovery/results/expanded_reference/common_gene_membership.csv'] : [])] : [])])])
const sourceFiles = computed(() => props.data.sources.filter(s => sourcePaths.value.includes(s.path)))
const coverageRows = computed(() => props.data.kinetics_summary.filter(r => r.scope === 'expanded' && r.comparison === 'available' &&
  r.dataset === dataset.value && r.branch === branch.value && r.parameter === parameter.value && rows.value.some(item => item.method === r.method)))
const expandedView = computed(() => mode.value === 'kinetics' && scope.value === 'expanded')
const format = value => value === null || value === undefined ? '—' : value.toFixed(4)
const bar = value => ({ left: `${value < 0 ? (value + 1) * 50 : 50}%`, width: `${Math.abs(value) * 50}%` })

</script>

<template>
  <div class="page-heading"><h1>Results</h1></div>
  <div class="tabs" role="tablist" aria-label="Evaluation level">
    <button v-for="[key, label] in [['cell','Cell-level'],['gene','Gene-level'],['kinetics','Kinetic-parameter recovery']]" :key="key"
      role="tab" :aria-selected="mode === key" :class="{ active: mode === key }" @click="mode = key">{{ label }}</button>
  </div>
  <div class="explorer">
    <aside class="filters panel">
      <div class="filter-title">Filters</div>
      <label for="dataset">Dataset</label><select id="dataset" v-model="dataset"><option v-for="d in datasets" :key="d.dataset_id" :value="d.dataset_id">{{ d.paper_name }}</option></select>
      <template v-if="mode === 'kinetics'">
        <label for="scope">Analysis scope</label><select id="scope" v-model="scope"><option value="expanded">Expanded reference genes</option><option value="full">Canonical scores</option><option value="selected">Figure 6 selected genes</option></select>
        <template v-if="scope === 'expanded'"><label for="comparison">Gene comparison</label><select id="comparison" v-model="comparison"><option v-for="(label,key) in comparisonNames" :key="key" :value="key">{{ label }}</option></select></template>
        <label v-if="dataset === 'sceu_organoid'" for="branch">Branch</label><select v-if="dataset === 'sceu_organoid'" id="branch" v-model="branch"><option v-for="b in branches" :key="b" :value="b">Branch {{ b }}</option></select>
        <label for="parameter">Kinetic parameter</label><select id="parameter" v-model="parameter"><option v-for="p in parameters" :key="p" :value="p">{{ p === 'alpha' ? 'α · transcription' : 'γ · degradation' }}</option></select>
      </template>
      <template v-if="mode !== 'cell'"><label for="gene">Gene / summary</label><select id="gene" v-model="gene">
        <option v-if="(mode === 'kinetics' && scope !== 'selected') || (mode === 'gene' && metric === 'cbdir')" value="__summary__">Mean</option>
        <option v-for="g in genes" :key="g" :value="g">{{ g }}</option></select></template>
      <label for="metric">Metric</label><select id="metric" v-model="metric"><option v-for="m in metrics" :key="m[0]" :value="m[0]">{{ m[1] }}</option></select>
      <label for="method-search">Search methods</label><input id="method-search" v-model="query" placeholder="Name or configuration ID">
      <details class="method-select"><summary>Methods <span>{{ selected.filter(m => availableMethods.includes(m)).length }} / {{ availableMethods.length }}</span></summary>
        <div class="method-actions"><button @click="selected = [...availableMethods]">Select all</button><button @click="selected = []">Clear</button></div>
        <label v-for="m in availableMethods" :key="m" class="checkbox"><input v-model="selected" type="checkbox" :value="m">{{ names[m] }}</label>
      </details>

    </aside>
    <section class="panel result-panel">
      <div class="result-heading"><div><p class="eyebrow">{{ datasetName }}</p><h2>{{ metricName }}</h2><p class="scope-text">{{ scopeLabel }}</p></div>
        </div>
      <div class="result-meta"><span>{{ finiteCount }} valid / {{ rows.length }} displayed methods</span><span v-if="mode !== 'cell'">{{ gene === '__summary__' ? '' : gene }}</span></div>
      <div class="table-toolbar"><div class="view-switch"><button class="text-button" :aria-pressed="display === 'table'" @click="display = 'table'">Table</button><button class="text-button" :aria-pressed="display === 'plot'" @click="display = 'plot'">Score plot</button><button v-if="mode !== 'cell' && gene !== '__summary__'" class="text-button" :aria-pressed="display === 'visual'" @click="display = 'visual'">{{ mode === 'gene' ? 'Phase portraits' : 'Trajectories' }}</button></div><label class="order-label">Order<select v-model="order" aria-label="Result order"><option value="method">Method</option><option value="descending">Score ↓</option><option value="ascending">Score ↑</option></select></label></div>
      <template v-if="display === 'visual' && rows.length"><GeneVisual v-if="visualIndex" :index="visualIndex" :colors="data.kinetics_colors" :level="mode" :dataset="dataset" :gene="gene" :branch="branch" :parameter="parameter" :rows="rows"/><p v-else-if="visualLoading" role="status">Loading visual index…</p><p v-else-if="visualError" role="alert">{{ visualError }}</p><p v-else class="empty">Visual assets unavailable. See Downloads for the source tables.</p></template>
      <div v-if="rows.length && display === 'plot'" class="comparison" aria-label="Method score comparison">
        <div class="axis"><span>−1</span><span>0</span><span>1</span></div>
        <div v-for="row in rows" :key="row.method" class="chart-row">
          <span class="chart-label" :title="row.name">{{ row.name }}</span>
          <div class="track"><div class="zero"></div><div v-if="row.value !== null" class="score-bar" :class="{ negative: row.value < 0 }" :style="bar(row.value)" :title="String(row.value)"></div></div>
          <span class="score-number">{{ format(row.value) }}</span>
        </div>
      </div>
      <div v-if="!rows.length" class="empty">No matching methods.</div>
      <div v-if="display === 'table'" class="table-scroll"><table><thead><tr><th>Method</th><th class="numeric">{{ metricName }}</th><th class="numeric">{{ mode === 'cell' ? 'Aligned cells' : gene === '__summary__' ? 'Valid genes' : mode === 'gene' ? 'CBDir evaluable cells' : 'Paired samples' }}</th><th>Scope / status</th></tr></thead>
        <tbody><tr v-for="row in rows" :key="row.method"><td>{{ row.name }}<small class="config-id">{{ row.method === 'dynamo_labelled' ? 'Labelled-data control' : row.method }}</small></td>
          <td class="numeric" :title="String(row.value)">{{ format(row.value) }}</td><td class="numeric">{{ row.n_genes ?? row.n_cells ?? row.n_samples ?? '—' }}</td>
          <td class="muted">{{ row.status || (mode === 'cell' ? 'Dataset score' : gene === '__summary__' ? 'Finite-gene mean' : 'Single gene') }}</td></tr></tbody></table></div>


      <details class="source-details"><summary>Source <span>{{ props.data.source_commit.slice(0,7) }}</span></summary>
        <div v-for="s in sourceFiles" :key="s.path" class="source-item"><a v-if="s.download" :href="`./${s.download}`" download>{{ s.path }} ↗</a><span v-else>{{ s.path }}</span><p>SHA-256 <code>{{ s.sha256 }}</code></p></div>
        <a href="./data/manifest.json" target="_blank" rel="noopener">Source manifest ↗</a>
      </details>
      <details v-if="expandedView" class="source-details"><summary>Coverage</summary>
        <p v-if="comparison !== 'available'" class="small muted">Common-gene membership is fixed for the comparison group. Coverage below refers to all available genes.</p>
        <div class="table-scroll"><table><thead><tr><th>Method</th><th class="numeric">Candidate</th><th class="numeric">Present</th><th class="numeric">Valid</th><th class="numeric">Undefined</th><th class="numeric">Unavailable</th></tr></thead><tbody>
          <tr v-for="row in coverageRows" :key="row.method"><td>{{ names[row.method] }}</td><td class="numeric">{{ row.candidate_n }}</td><td class="numeric">{{ row.present_n }}</td><td class="numeric">{{ row.n_genes }}</td><td class="numeric">{{ row.undefined_n }}</td><td class="numeric">{{ row.unavailable_n }}</td></tr>
        </tbody></table></div>
      </details>
    </section>
  </div>
</template>
