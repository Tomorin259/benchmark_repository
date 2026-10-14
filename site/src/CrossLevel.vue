<script setup>
import { computed, ref } from 'vue'
import ImageViewer from './ImageViewer.vue'
const props = defineProps({ data:Object, caseId:String })
const emit = defineEmits(['update:caseId'])
const enlarged = ref(null)
const current = computed(()=>props.data.cross_cases.find(c=>c.id===props.caseId)||props.data.cross_cases[0])
const names = Object.fromEntries(props.data.methods.map(m=>[m.config_id,m.manuscript_name]))
const rows = computed(()=>props.data.methods.map(m=>props.data.cross_level.find(r=>r.case===current.value.id&&r.method===m.config_id)).filter(Boolean))
const dataset = computed(()=>props.data.datasets.find(d=>d.dataset_id===current.value.dataset))
const format = value=>value==null?'—':Number(value).toFixed(4)
</script>

<template>
  <div class="explorer"><aside class="filters panel"><label for="cross-case">Case</label><select id="cross-case" :value="current.id" @change="emit('update:caseId',$event.target.value)"><option v-for="c in data.cross_cases" :key="c.id" :value="c.id">{{c.label}}</option></select></aside>
    <section class="panel result-panel"><div class="result-heading"><h2>{{current.label}}</h2><span class="tag">{{current.n_genes}} genes</span></div>
      <div class="table-scroll cross-table"><table><thead><tr><th rowspan="2">Method</th><th colspan="2" class="numeric">Cell-level</th><th colspan="2" class="numeric">Gene-level · Mean ± SD</th></tr><tr><th class="numeric">CBDir</th><th class="numeric">BATC</th><th class="numeric">CBDir-gene</th><th class="numeric">CosSim</th></tr></thead><tbody><tr v-for="r in rows" :key="r.method"><td>{{names[r.method]}}</td><td class="numeric">{{format(r.cell_cbdir)}}</td><td class="numeric">{{format(r.cell_batc)}}</td><td class="numeric">{{format(r.gene_cbdir_mean)}}<small class="summary-sd">± {{format(r.gene_cbdir_sd)}}</small></td><td class="numeric">{{format(r.gene_cossim_mean)}}<small class="summary-sd">± {{format(r.gene_cossim_sd)}}</small></td></tr></tbody></table></div>
      <p class="small muted compact-note">SD across evaluated genes.</p>
      <div class="guide-links"><a :href="`#results?level=cell&dataset=${current.dataset}`">{{dataset.paper_name}} · Cell →</a><a :href="`#results?level=gene&dataset=${current.dataset}`">Gene →</a></div>
      <details class="source-details"><summary>All simulated comparisons</summary><button class="phase-image" aria-label="Enlarge all simulated cross-level comparisons" @click="enlarged={...data.result_images['simulation-cross-level'],title:'All simulated comparisons'}"><img :src="'./'+data.result_images['simulation-cross-level'].path" alt="All four simulated Cell/Gene comparisons" loading="lazy"></button></details>
      <p class="reference-link"><a href="#downloads">Download paired results →</a></p>
    </section>
  </div>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>
