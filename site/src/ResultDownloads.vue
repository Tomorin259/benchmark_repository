<script setup>
import { computed, ref, watch } from 'vue'

const props = defineProps({ data: { type: Object, required: true } })
const dataset = ref('all'), method = ref('all'), kind = ref('all'), page = ref(1)
const kinds = [
  ['method_outputs', 'Native method outputs'], ['kinetics_outputs', 'Kinetics outputs'],
  ['simulation_outputs', 'Simulation outputs'], ['ensemble_outputs', 'Ensemble outputs'],
  ['diagnostics', 'TIvelo / MultiVelo controls'], ['evaluation_assets', 'Reference moments'],
]
const packages = computed(() => props.data.result_packages.map(item => ({
  ...item,
  kind: item.repository_path.split('/')[1],
  datasetIds: item.dataset_ids.split(';').filter(Boolean),
  methodIds: item.config_ids.split(';').filter(Boolean),
  datasetLabel: item.dataset_ids.split(';').map(id => props.data.datasets.find(d => d.dataset_id === id)?.paper_name || id).join(', '),
  methodLabel: item.config_ids.split(';').filter(Boolean).map(id => props.data.methods.find(m => m.config_id === id)?.manuscript_name || id).join(', '),
})))
const filtered = computed(() => packages.value.filter(item =>
  (dataset.value === 'all' || item.datasetIds.includes(dataset.value)) &&
  (method.value === 'all' || item.methodIds.includes(method.value)) &&
  (kind.value === 'all' || item.kind === kind.value)))
const pageCount = computed(() => Math.max(1, Math.ceil(filtered.value.length / 25)))
const visible = computed(() => filtered.value.slice((page.value - 1) * 25, page.value * 25))
const totalBytes = computed(() => packages.value.reduce((sum, item) => sum + item.bytes, 0))
const size = bytes => bytes >= 1e9 ? `${(bytes / 1e9).toFixed(2)} GB` : bytes >= 1e6 ? `${(bytes / 1e6).toFixed(1)} MB` : `${(bytes / 1e3).toFixed(1)} KB`
const kindName = id => kinds.find(([key]) => key === id)?.[1] || id
watch([dataset, method, kind], () => { page.value = 1 })
</script>

<template>
  <section class="panel downloads-section" aria-labelledby="result-download-heading">
    <h2 id="result-download-heading">Full method outputs</h2>
    <p class="small muted">{{packages[0].release_version}} · {{packages.length}} packages · {{size(totalBytes)}}. Native H5AD outputs are grouped in lossless TAR.GZ archives; reference moments are provided as NPZ files.</p>
    <div class="guide-links">
      <a href="https://huggingface.co/datasets/Zelin259/rna-velocity-benchmark/tree/v1.0.0" target="_blank" rel="noopener">Hugging Face release ↗</a>
      <a href="./data/files/data/RESULT_DOWNLOADS.md" download>Extraction guide</a>
      <a href="./data/files/data/result_package_manifest.tsv" download>Package manifest</a>
      <a href="./data/files/data/result_member_manifest.tsv" download>Original file manifest</a>
    </div>
    <div class="classification-filters">
      <label>Dataset<select v-model="dataset" aria-label="Result download dataset"><option value="all">All datasets</option><option v-for="d in data.datasets" :key="d.dataset_id" :value="d.dataset_id">{{d.paper_name}}</option></select></label>
      <label>Method<select v-model="method" aria-label="Result download method"><option value="all">All methods</option><option v-for="m in data.methods" :key="m.config_id" :value="m.config_id">{{m.manuscript_name}}</option></select></label>
      <label>Output type<select v-model="kind" aria-label="Result download type"><option value="all">All output types</option><option v-for="[id,label] in kinds" :key="id" :value="id">{{label}}</option></select></label>
    </div>
    <p class="small muted" aria-live="polite">{{filtered.length}} matching {{filtered.length === 1 ? 'package' : 'packages'}}. Download the dataset and method outputs you need.</p>
    <div v-for="file in visible" :key="file.package_id" class="download-row result-package-row">
      <div>
        <strong>{{file.datasetLabel}} · {{file.methodLabel || kindName(file.kind)}}</strong>
        <span class="small muted">{{kindName(file.kind)}} · {{file.format.toUpperCase()}} · {{size(file.bytes)}} · {{file.scientific_files}} original {{file.scientific_files === 1 ? 'file' : 'files'}}</span>
        <span class="small package-path">{{file.repository_path.split('/').slice(1).join('/')}}</span>
      </div>
      <a class="button secondary" :href="file.url" target="_blank" rel="noopener">Download {{file.format.toUpperCase()}}</a>
    </div>
    <p v-if="!filtered.length" class="empty">No packages match these filters.</p>
    <div v-if="filtered.length > 25" class="pagination">
      <span>Page {{page}} of {{pageCount}}</span>
      <button class="button secondary" :disabled="page === 1" @click="page--">Previous</button>
      <button class="button secondary" :disabled="page === pageCount" @click="page++">Next</button>
    </div>
    <p class="small muted">Extract each TAR.GZ package into a separate folder. NPZ files can be used directly.</p>
  </section>
</template>

<style scoped>
.result-package-row>div{min-width:0}.package-path{overflow-wrap:anywhere}.result-package-row>a{flex-shrink:0}.pagination{margin-bottom:16px}
@media(max-width:760px){.result-package-row{flex-wrap:wrap}.result-package-row>div{flex-basis:100%}}
</style>
