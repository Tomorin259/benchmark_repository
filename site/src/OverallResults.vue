<script setup>
import { computed, ref } from 'vue'
import ImageViewer from './ImageViewer.vue'

const props = defineProps({ data:Object, cellOnly:Boolean, task:{type:String,default:'cell'}, statistic:{type:String,default:'mean'}, view:{type:String,default:'table'} })
const emit = defineEmits(['update:task','update:statistic','update:view'])
const opened = ref(''), enlarged = ref(null)
const names = Object.fromEntries(props.data.methods.map(m=>[m.config_id,m.manuscript_name]))
names.dynamo_labelled = 'Dynamo (scEU-seq)'
const cellMetrics = [['consistency','Consistency'],['iccoh','ICCoh'],['cbdir','CBDir'],['transcosine','TransCosine'],['batc','BATC']]
const metrics = computed(()=>props.cellOnly||props.task==='cell'?cellMetrics:props.task==='gene'?
  [['Real-CBDir','Real · CBDir-gene'],['Sim-CBDir','Simulation · CBDir-gene'],['Sim-CosSim','Simulation · CosSim']]:
  [['Kinetic-Spearman-alpha','α · Spearman'],['Kinetic-Spearman-gamma','γ · Spearman']])
const cellTask = computed(()=>props.cellOnly||props.task==='cell')
const sortMetric = ref(''), descending = ref(true)
const format = value=>value==null?'—':Number(value).toFixed(4)
const rows = computed(()=>{
  const source = cellTask.value?props.data.cell_summary:props.data.overall_summary
  const result = [...props.data.methods.map(m=>m.config_id),'dynamo_labelled'].map(method=>({
    method, values: Object.fromEntries(source.filter(r=>r.method===method&&metrics.value.some(m=>m[0]===r.metric)).map(r=>[r.metric,r])),
  })).filter(r=>Object.keys(r.values).length)
  if(sortMetric.value&&metrics.value.some(m=>m[0]===sortMetric.value))result.sort((a,b)=>{
    const field=cellTask.value?props.statistic:'mean', av=a.values[sortMetric.value]?.[field], bv=b.values[sortMetric.value]?.[field]
    if(av==null)return bv==null?0:1
    if(bv==null)return -1
    return (descending.value?-1:1)*(av-bv)
  })
  return result
})
const datasets = computed(()=>props.data.datasets.filter(d=>props.data.cell.some(r=>r.dataset===d.dataset_id)))
function sort(key){if(sortMetric.value===key)descending.value=!descending.value;else{sortMetric.value=key;descending.value=true}}
function changeTask(task){opened.value='';sortMetric.value='';emit('update:task',task);emit('update:view','table')}
const detailRows = id=>datasets.value.map(d=>({dataset:d,...props.data.cell.find(r=>r.dataset===d.dataset_id&&r.method===id)}))
</script>

<template>
  <section class="panel result-panel overall-panel">
    <div class="result-heading"><h2>{{cellOnly?'All real datasets':'Benchmark summary'}}</h2><span v-if="cellTask" class="tag">14 datasets · Equal weight</span></div>
    <div class="table-toolbar">
      <div v-if="!cellOnly" class="view-switch" aria-label="Summary task"><button v-for="[key,label] in [['cell','Cell'],['gene','Gene'],['kinetics','Kinetics']]" :key="key" class="text-button" :aria-pressed="task===key" @click="changeTask(key)">{{label}}</button></div>
      <div v-if="cellTask" class="table-controls"><label class="order-label" for="summary-statistic">Statistic<select id="summary-statistic" :value="statistic" @change="emit('update:statistic',$event.target.value)"><option value="mean">Mean</option><option value="median">Median</option></select></label><div class="view-switch"><button class="text-button" :aria-pressed="view!=='distribution'" @click="emit('update:view','table')">Table</button><button class="text-button" :aria-pressed="view==='distribution'" @click="emit('update:view','distribution')">Distribution</button></div></div>
    </div>
    <p v-if="!cellTask&&task==='gene'" class="small muted compact-note">Each method’s available finite gene scores are averaged within each dataset or scenario, then the seven real datasets or three simulation scenarios receive equal weight. Values show mean ± SD across these means. <a href="#results?level=gene">Per-dataset scores and gene counts →</a></p>
    <p v-if="!cellTask&&task==='kinetics'" class="small muted compact-note">Finite gene scores are averaged within each branch. Organoid branches and available systems receive equal weight. Values show mean ± SD across system means; SD requires at least two systems.</p>
    <div v-if="!cellTask||view!=='distribution'" class="table-scroll"><table class="summary-table"><thead><tr><th>Method</th><th v-for="[key,label] in metrics" :key="key" class="numeric"><button class="sort-heading" @click="sort(key)">{{label}}{{sortMetric===key?(descending?' ↓':' ↑'):''}}</button></th></tr></thead>
      <tbody><template v-for="row in rows" :key="row.method"><tr :class="{selected:opened===row.method}"><td><button v-if="cellTask" class="method-link" :aria-expanded="opened===row.method" @click="opened=opened===row.method?'':row.method">{{names[row.method]}}</button><span v-else>{{names[row.method]}}</span></td><td v-for="[key] in metrics" :key="key" class="numeric" :title="row.values[key]?'SD '+format(row.values[key].sample_sd):undefined">{{format(row.values[key]?.[cellTask?statistic:'mean'])}}<small v-if="!cellTask&&row.values[key]?.sample_sd!=null" class="summary-sd">± {{format(row.values[key].sample_sd)}}</small><small v-if="!cellTask&&task==='kinetics'&&row.values[key]" class="summary-sd">{{row.values[key].n_components}} {{row.values[key].n_components===1?'system':'systems'}}</small></td></tr>
        <tr v-if="cellTask&&opened===row.method"><td :colspan="metrics.length+1" class="method-detail"><div class="table-scroll"><table><thead><tr><th>Dataset</th><th v-for="[key,label] in cellMetrics" :key="key" class="numeric">{{label}}</th></tr></thead><tbody><tr v-for="detail in detailRows(row.method)" :key="detail.dataset.dataset_id"><td><a :href="`#results?level=cell&dataset=${detail.dataset.dataset_id}&method=${row.method}`">{{detail.dataset.paper_name}} →</a></td><td v-for="[key] in cellMetrics" :key="key" class="numeric">{{format(detail[key])}}</td></tr></tbody></table></div></td></tr>
      </template></tbody></table></div>
    <div v-else class="distribution-grid"><article v-for="[key,label] in cellMetrics" :key="key"><h3>{{label}}</h3><button class="phase-image" :aria-label="`Enlarge ${label} distribution`" @click="enlarged={...data.result_images['cell-'+key],title:label+' · 14 datasets'}"><img :src="'./'+data.result_images['cell-'+key].path" :alt="label+' across 14 real datasets'" loading="lazy"></button></article></div>
  <p class="reference-link"><a href="#downloads">Source data →</a></p>
    <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
  </section>
</template>
