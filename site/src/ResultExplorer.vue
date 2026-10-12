<script setup>
import { computed, nextTick, onMounted, ref, watch } from 'vue'
import GeneVisual from './GeneVisual.vue'

const props = defineProps({ data:Object, initialLevel:{type:String,default:'cell'} })
const mode=ref(['cell','gene','kinetics'].includes(props.initialLevel)?props.initialLevel:'cell')
const dataset=ref(mode.value==='kinetics'?'sceu_rpe1':'pancreas'), metric=ref(mode.value==='kinetics'?'score':'cbdir')
const scope=ref('full'), comparison=ref('available'), branch=ref(''), parameter=ref('alpha'), gene=ref('__summary__')
const query=ref(''), order=ref('method'), display=ref('table'), selected=ref(props.data.methods.map(m=>m.config_id).concat('dynamo_labelled'))
const visualIndex=ref(null), visualLoading=ref(false), visualError=ref(''), ready=ref(false)
const streamMethods=ref([]), compareOpen=ref(false), enlarged=ref(null), dialog=ref(null), copied=ref(false)
const names=Object.fromEntries(props.data.methods.map(m=>[m.config_id,m.manuscript_name]))
names.dynamo_labelled='Dynamo (scEU-seq)'
const cellMetrics=[['consistency','Consistency'],['iccoh','ICCoh'],['cbdir','CBDir'],['transcosine','TransCosine'],['batc','BATC']]
const comparisonNames={available:'Each method’s available genes',core_conventional:'Fixed common set · core RNA methods',all_conventional:'Fixed common set · all RNA methods',all_displayed_methods:'Fixed common set · RNA methods + scEU comparator'}
const definitions={consistency:'Local agreement between neighboring velocity vectors.',iccoh:'Velocity coherence within annotated cell populations.',cbdir:'Directional agreement across annotated transitions between cell states.',transcosine:'Agreement between the velocity transition graph and annotated transitions.',batc:'Directional agreement along the broader reference lineage.',cossim:'Agreement between inferred and ground-truth gene velocities in simulations.',score:'Agreement between predicted and experimentally fitted kinetic profiles across paired trajectory positions.'}
const tables=computed(()=>props.data[mode.value])
const datasets=computed(()=>props.data.datasets.filter(d=>tables.value.some(r=>r.dataset===d.dataset_id)))
const datasetName=computed(()=>props.data.datasets.find(d=>d.dataset_id===dataset.value)?.paper_name)
const branches=computed(()=>[...new Set(props.data.kinetics.filter(r=>r.dataset===dataset.value&&r.scope===scope.value).map(r=>r.branch))].sort())
const parameters=computed(()=>[...new Set(props.data.kinetics.filter(r=>r.dataset===dataset.value&&r.scope===scope.value&&r.branch===branch.value).map(r=>r.parameter))].sort())
const commonGenes=computed(()=>props.data.kinetics_common_membership.find(r=>r.dataset===dataset.value&&r.branch===branch.value&&r.parameter===parameter.value&&r.comparison===comparison.value)?.genes??[])
const commonMethods=computed(()=>props.data.kinetics_summary.filter(r=>r.scope==='expanded'&&r.dataset===dataset.value&&r.branch===branch.value&&r.parameter===parameter.value&&r.comparison===comparison.value).map(r=>r.method))
const candidates=computed(()=>tables.value.filter(r=>r.dataset===dataset.value&&(mode.value!=='kinetics'||(r.scope===scope.value&&r.branch===branch.value&&r.parameter===parameter.value&&(scope.value!=='expanded'||comparison.value==='available'||commonGenes.value.includes(r.gene)&&commonMethods.value.includes(r.method))))))
const genes=computed(()=>[...new Set(candidates.value.map(r=>r.gene).filter(Boolean))].sort())
const metrics=computed(()=>mode.value==='cell'?cellMetrics:mode.value==='gene'?[['cbdir','CBDir-gene'],...(dataset.value.startsWith('sim_')?[['cossim','CosSim']]:[])]:[['score','Spearman correlation']])
const metricName=computed(()=>metrics.value.find(m=>m[0]===metric.value)?.[1])
const summaryAllowed=computed(()=>mode.value==='kinetics'?scope.value!=='selected':metric.value==='cbdir')
watch(mode,value=>{
  dataset.value=value==='kinetics'?'sceu_rpe1':'pancreas'; metric.value=value==='kinetics'?'score':'cbdir'
  query.value=''; gene.value='__summary__'; display.value='table'; streamMethods.value=[]
})
watch(()=>props.initialLevel,value=>{if(['cell','gene','kinetics'].includes(value)&&value!==mode.value)mode.value=value})
watch(branches,value=>{if(!value.includes(branch.value))branch.value=value[0]??''},{immediate:true})
watch(parameters,value=>{if(!value.includes(parameter.value))parameter.value=value[0]??'alpha'},{immediate:true})
watch(scope,value=>{if(value!=='expanded')comparison.value='available'})
watch([genes,summaryAllowed],()=>{
  if(!genes.value.includes(gene.value)&&!(gene.value==='__summary__'&&summaryAllowed.value))gene.value=gene.value==='__summary__'&&summaryAllowed.value?'__summary__':genes.value[0]??''
},{immediate:true})
watch(metrics,value=>{if(!value.some(m=>m[0]===metric.value))metric.value=value[0]?.[0]},{immediate:true})
watch(gene,value=>{if(value==='__summary__')display.value='table'})
const unfilteredRows=computed(()=>{
  const records=mode.value==='kinetics'&&gene.value==='__summary__'?props.data.kinetics_summary.filter(r=>r.dataset===dataset.value&&r.branch===branch.value&&r.parameter===parameter.value&&r.scope===scope.value&&r.comparison===comparison.value)
    :mode.value==='gene'&&gene.value==='__summary__'?props.data.gene_summary.filter(r=>r.dataset===dataset.value)
    :candidates.value.filter(r=>mode.value==='cell'||r.gene===gene.value)
  return records.map(r=>({...r,name:names[r.method],value:r[metric.value]}))
})
const availableMethods=computed(()=>[...new Set(unfilteredRows.value.map(r=>r.method))])
const rows=computed(()=>unfilteredRows.value.filter(r=>mode.value==='cell'||selected.value.includes(r.method)&&r.name.toLowerCase().includes(query.value.toLowerCase())).sort((a,b)=>{
  const ids=mode.value==='kinetics'?props.data.kinetics_method_order:props.data.methods.map(m=>m.config_id)
  if(order.value==='method')return ids.indexOf(a.method)-ids.indexOf(b.method)
  if(a.value===null)return b.value===null?ids.indexOf(a.method)-ids.indexOf(b.method):1
  if(b.value===null)return -1
  return (order.value==='ascending'?a.value-b.value:b.value-a.value)||ids.indexOf(a.method)-ids.indexOf(b.method)
}))
const scopeLabel=computed(()=>mode.value==='cell'?props.data.dataset_descriptions[dataset.value]:mode.value==='gene'?metric.value==='cossim'?'Simulation ground-truth comparison':gene.value==='__summary__'?'Mean across evaluated genes with finite CBDir scores':'Single-gene CBDir on the fixed reference neighborhood':scope.value==='selected'?`Figure 6 example genes · ${dataset.value==='sceu_rpe1'?19:branch.value==='1'?6:8} genes`:scope.value==='full'?'Main benchmark gene scores · source for Figure 7 summaries':`Extended reference comparison · ${props.data.kinetics_expanded_scope[dataset.value]} candidate genes`)
const coverageRows=computed(()=>props.data.kinetics_summary.filter(r=>r.scope==='expanded'&&r.comparison==='available'&&r.dataset===dataset.value&&r.branch===branch.value&&r.parameter===parameter.value&&rows.value.some(x=>x.method===r.method)))
const streamRows=computed(()=>streamMethods.value.map(id=>rows.value.find(r=>r.method===id)).filter(Boolean).map(r=>({...r,asset:props.data.streams.find(s=>s.dataset===dataset.value&&s.method===r.method)})))
const previewChoices=computed(()=>rows.value.map(r=>({...r,asset:props.data.streams.find(s=>s.dataset===dataset.value&&s.method===r.method)})).filter(r=>r.asset))
const previewMethod=computed({get:()=>streamMethods.value[0]??'',set:value=>{if(value)streamMethods.value=[value]}})
watch([mode,rows],([level,value])=>{
  const valid=streamMethods.value.filter(id=>value.some(r=>r.method===id))
  streamMethods.value=valid.length?valid:(level==='cell'&&value.length?[value[0].method]:[])
},{immediate:true})
const format=value=>value===null||value===undefined?'—':value.toFixed(4)
const status=row=>metric.value==='cossim'?'':(row.status&&row.status!=='ok'?row.status.replaceAll('_',' '):'')
const hasStatus=computed(()=>rows.value.some(r=>status(r)))
const bar=value=>({left:`${value<0?(value+1)*50:50}%`,width:`${Math.abs(value)*50}%`})
const finiteCount=computed(()=>rows.value.filter(r=>r.value!==null).length)
function chooseSubject(value){gene.value=value==='summary'?'__summary__':genes.value[0]??'';display.value=value==='summary'?'table':'visual'}
function sortMetric(key){if(metric.value===key)order.value=order.value==='descending'?'ascending':'descending';else{metric.value=key;order.value='descending'}}
async function preview(id){streamMethods.value=[id];compareOpen.value=false;display.value='table';await nextTick();document.getElementById('stream-preview')?.scrollIntoView({block:'start',behavior:'smooth'})}
function toggleCompare(id){if(streamMethods.value.includes(id)){if(streamMethods.value.length>1)streamMethods.value=streamMethods.value.filter(m=>m!==id);return}if(streamMethods.value.length<4)streamMethods.value=[...streamMethods.value,id]}
function stepStream(direction){const i=rows.value.findIndex(r=>r.method===streamMethods.value[0]);const row=rows.value[(i+direction+rows.value.length)%rows.value.length];if(row)streamMethods.value=[row.method]}
watch(enlarged,async value=>{await nextTick();if(value)dialog.value?.showModal()})
watch([display,gene],async()=>{
  if(display.value!=='visual'||gene.value==='__summary__'||visualIndex.value||visualLoading.value||!props.data.visuals)return
  visualLoading.value=true;visualError.value=''
  try{const response=await fetch(`./${props.data.visuals.index}`);if(!response.ok)throw Error(response.status);visualIndex.value=await response.json()}
  catch(e){console.error(e);visualError.value='The visual could not be loaded. Choose Table or try another gene.'}
  finally{visualLoading.value=false}
})
const shareUrl=computed(()=>{
  const p=new URLSearchParams({level:mode.value,dataset:dataset.value,metric:metric.value,order:order.value,view:display.value})
  if(mode.value!=='cell')p.set('gene',gene.value)
  if(mode.value==='kinetics'){p.set('scope',scope.value);p.set('comparison',comparison.value);p.set('branch',branch.value);p.set('parameter',parameter.value)}
  if(mode.value!=='cell'&&query.value)p.set('search',query.value)
  if(mode.value!=='cell'&&selected.value.length<props.data.methods.length+1)p.set('methods',selected.value.join(','))
  if(mode.value==='cell'&&streamMethods.value.length===1)p.set('method',streamMethods.value[0])
  if(mode.value==='cell'&&streamMethods.value.length>1)p.set('compare',streamMethods.value.join(','))
  return `${location.origin}${location.pathname}#results?${p}`
})
watch(shareUrl,value=>{if(ready.value)history.replaceState(null,'',value)})
async function copyLink(){try{await navigator.clipboard.writeText(shareUrl.value);copied.value=true;setTimeout(()=>copied.value=false,2000)}catch{copied.value=false}}
onMounted(async()=>{
  const p=new URLSearchParams(location.hash.split('?')[1]);await nextTick()
  if(datasets.value.some(d=>d.dataset_id===p.get('dataset')))dataset.value=p.get('dataset')
  if(['full','selected','expanded'].includes(p.get('scope')))scope.value=p.get('scope')
  await nextTick();if(branches.value.includes(p.get('branch')))branch.value=p.get('branch')
  await nextTick();if(parameters.value.includes(p.get('parameter')))parameter.value=p.get('parameter')
  if(comparisonNames[p.get('comparison')])comparison.value=p.get('comparison')
  if(metrics.value.some(m=>m[0]===p.get('metric')))metric.value=p.get('metric')
  await nextTick();if(genes.value.includes(p.get('gene'))||p.get('gene')==='__summary__'&&summaryAllowed.value)gene.value=p.get('gene')
  if(['method','ascending','descending'].includes(p.get('order')))order.value=p.get('order')
  if(p.has('methods'))selected.value=p.get('methods').split(',').filter(m=>names[m])
  query.value=mode.value==='cell'?'':p.get('search')??''
  await nextTick();if(['table','plot','visual'].includes(p.get('view'))&&!(gene.value==='__summary__'&&p.get('view')==='visual'))display.value=p.get('view')
  if(mode.value==='cell'){
    const requested=p.has('compare')?p.get('compare').split(','):p.has('method')?[p.get('method')]:p.has('streams')?p.get('streams').split(','):[]
    const valid=requested.filter(m=>rows.value.some(r=>r.method===m)).slice(0,4)
    if(valid.length)streamMethods.value=valid
  }
  ready.value=true
  history.replaceState(null,'',shareUrl.value)
})
</script>

<template>
  <div class="page-heading"><h1>Results</h1><button class="button secondary" @click="copyLink">{{ copied?'Link copied':'Share this view' }}</button></div>
  <div class="tabs" aria-label="Evaluation tasks"><button v-for="[key,label] in [['cell','Cell-level'],['gene','Gene-level'],['kinetics','Kinetic-parameter recovery']]" :key="key" :class="{active:mode===key}" :aria-pressed="mode===key" @click="mode=key">{{ label }}</button></div>
  <div class="explorer"><aside class="filters panel"><div class="filter-title">Choose a comparison</div>
    <label for="dataset">{{ mode==='kinetics'?'System':'Dataset' }}</label><select id="dataset" v-model="dataset"><option v-for="d in datasets" :key="d.dataset_id" :value="d.dataset_id">{{ d.paper_name }}</option></select>
    <template v-if="mode==='kinetics'"><label for="scope">Analysis</label><select id="scope" v-model="scope"><option value="full">Main benchmark scores (Figure 7)</option><option value="selected">Figure 6 example genes</option><option value="expanded">Extended reference-gene analysis</option></select>
      <template v-if="scope==='expanded'"><label for="comparison">Gene comparison</label><select id="comparison" v-model="comparison"><option v-for="(label,key) in comparisonNames" :key="key" :value="key">{{ label }}</option></select></template>
      <template v-if="dataset==='sceu_organoid'"><label for="branch">Differentiation branch</label><select id="branch" v-model="branch"><option value="1">Branch 1 · Secretory</option><option value="2">Branch 2 · Enterocyte</option></select></template>
      <label for="parameter">Parameter</label><select id="parameter" v-model="parameter"><option v-for="p in parameters" :key="p" :value="p">{{ p==='alpha'?'α · Transcription':'γ · Degradation' }}</option></select>
    </template>
    <template v-if="mode!=='cell'"><label>Evaluation object</label><div class="subject-switch"><button :disabled="!summaryAllowed" :aria-pressed="gene==='__summary__'" @click="chooseSubject('summary')">Summary</button><button :aria-pressed="gene!=='__summary__'" @click="chooseSubject('single')">Single gene</button></div><template v-if="gene!=='__summary__'"><label for="gene">Gene</label><select id="gene" v-model="gene"><option v-for="g in genes" :key="g">{{ g }}</option></select></template></template>
    <template v-if="mode!=='kinetics'"><label for="metric">{{ mode==='cell'?'Sort / plot metric':'Metric' }}</label><select id="metric" v-model="metric"><option v-for="m in metrics" :key="m[0]" :value="m[0]">{{ m[1] }}</option></select></template>
    <template v-if="mode!=='cell'"><label for="method-search">Search methods</label><input id="method-search" v-model="query" placeholder="Method name">
      <details class="method-select"><summary>Methods <span>{{ rows.length }} / {{ availableMethods.length }}</span></summary><div class="method-actions"><button @click="selected=[...availableMethods]">Select all</button><button @click="selected=[]">Clear</button></div><label v-for="m in availableMethods" :key="m" class="checkbox"><input v-model="selected" type="checkbox" :value="m">{{ names[m] }}</label></details>
    </template>
  </aside>
  <section class="panel result-panel"><div class="result-heading"><div><p class="eyebrow">{{ datasetName }}{{ branch?` · ${branch==='1'?'Secretory':'Enterocyte'} branch`:'' }}</p><h2>{{ mode==='cell'?'Method scores':gene==='__summary__'?'Gene-score summary':gene }}</h2><p class="scope-text">{{ scopeLabel }}</p></div><a class="text-link" href="#directory">Data & methods →</a></div>
    <p class="metric-description"><strong>{{ metricName }}:</strong> {{ mode==='gene'&&metric==='cbdir'?'Gene-specific directional agreement across annotated cell-state transitions.':definitions[metric] }} Higher values indicate stronger agreement.</p>
    <p v-if="mode==='kinetics'&&scope==='expanded'&&comparison!=='available'" class="small notice">The common-gene set is frozen for {{ commonMethods.map(m=>names[m]).join(', ') }}. Filtering methods changes the display, not the gene set.</p>
    <section v-if="mode==='cell'&&streamRows.length" id="stream-preview" class="stream-preview"><div class="stream-toolbar"><div><p class="eyebrow">Cell-level preview</p><h2>{{ datasetName }} velocity field</h2></div><label class="preview-method">Method<select v-model="previewMethod"><option v-for="r in previewChoices" :key="r.method" :value="r.method">{{ r.name }}</option></select></label><div><button class="text-button" v-if="streamRows.length===1" @click="stepStream(-1)">Previous</button><button class="text-button" v-if="streamRows.length===1" @click="stepStream(1)">Next</button><button class="text-button" @click="compareOpen=!compareOpen">{{ compareOpen?'Hide method list':'Compare methods' }}</button></div></div>
      <fieldset v-if="compareOpen" class="compare-picker"><legend>Methods to show (up to four)</legend><label v-for="r in previewChoices" :key="r.method"><input type="checkbox" :checked="streamMethods.includes(r.method)" :disabled="streamMethods.length>=4&&!streamMethods.includes(r.method)" @change="toggleCompare(r.method)">{{ r.name }}</label></fieldset>
      <div class="stream-grid" :class="{single:streamRows.length===1}"><article v-for="r in streamRows" :key="r.method"><h3>{{ r.name }}</h3><p class="small">{{ metricName }} {{ format(r.value) }}</p><button v-if="r.asset" class="phase-image" :aria-label="`Enlarge ${r.name} streamlines`" @click="enlarged=r"><img :src="`./${r.asset.path}`" :alt="`${datasetName}: ${r.name} cell velocity field`" loading="lazy"></button><p v-else>Streamlines unavailable.</p><dl v-if="streamRows.length===1" class="metric-values"><div v-for="[key,label] in cellMetrics" :key="key"><dt>{{ label }}</dt><dd>{{ format(r[key]) }}</dd></div></dl><details v-if="r.n_cells<r.n_reference" class="coverage-note"><summary>Partial cell coverage</summary><p>{{ r.n_cells }} of {{ r.n_reference }} reference cells are aligned for this method. This is cell coverage, not the valid-sample count of each metric.</p></details></article></div>
      <div class="cell-legend" v-if="streamRows[0]?.asset"><span v-for="item in streamRows[0].asset.legend" :key="item.label"><i :style="{background:item.color}"></i>{{ item.label }}</span></div><p class="small muted">The field is shown on the dataset reference coordinates. Use the scores below to compare agreement with the evaluation references.</p>
    </section>
    <div class="table-toolbar"><div class="view-switch"><button class="text-button" :aria-pressed="display==='table'" @click="display='table'">{{ mode==='cell'?'All metrics':'Table' }}</button><button class="text-button" :aria-pressed="display==='plot'" @click="display='plot'">Score plot</button><button v-if="mode!=='cell'" class="text-button" :disabled="gene==='__summary__'" :aria-pressed="display==='visual'" @click="display='visual'">{{ mode==='gene'?'Phase portraits':'Trajectories' }}</button></div><label class="order-label">Order<select v-model="order"><option value="method">Method</option><option value="descending">Score ↓</option><option value="ascending">Score ↑</option></select></label></div>
    <p v-if="mode!=='cell'&&gene==='__summary__'" class="small muted">Choose Single gene to inspect {{ mode==='gene'?'phase portraits':'reference and predicted trajectories' }}.</p>
    <div v-if="display==='table'&&mode==='cell'" class="table-scroll"><table class="cell-scores"><thead><tr><th>Method</th><th class="numeric" v-for="[key,label] in cellMetrics" :key="key"><button @click="sortMetric(key)" :aria-label="`Sort by ${label}`">{{ label }}{{ metric===key&&order!=='method'?(order==='descending'?' ↓':' ↑'):'' }}</button></th></tr></thead><tbody><tr v-for="r in rows" :key="r.method"><td><button class="method-link" @click="preview(r.method)">{{ r.name }}</button></td><td class="numeric" v-for="[key] in cellMetrics" :key="key">{{ format(r[key]) }}</td></tr></tbody></table></div>
    <div v-else-if="display==='table'" class="table-scroll"><table><thead><tr><th>Method</th><th class="numeric">{{ metricName }}</th><th v-if="metric!=='cossim'" class="numeric">{{ gene==='__summary__'?'Valid genes':mode==='gene'?'Evaluable cells':'Paired positions' }}</th><th v-if="hasStatus">Status</th></tr></thead><tbody><tr v-for="r in rows" :key="r.method"><td>{{ r.name }}<small v-if="r.method==='dynamo_labelled'" class="muted">Labeling-data comparator</small></td><td class="numeric" :title="String(r.value)">{{ format(r.value) }}</td><td v-if="metric!=='cossim'" class="numeric">{{ r.n_genes??r.n_cells??r.n_samples??'—' }}</td><td v-if="hasStatus">{{ status(r) }}</td></tr></tbody></table></div>
    <div v-if="display==='plot'&&rows.length" class="comparison"><div class="axis"><span>−1</span><span>0</span><span>1</span></div><div class="chart-row" v-for="r in rows" :key="r.method"><span class="chart-label">{{ r.name }}</span><div class="track"><div class="zero"></div><div v-if="r.value!==null" class="score-bar" :class="{negative:r.value<0}" :style="bar(r.value)"></div></div><span class="score-number">{{ format(r.value) }}</span></div></div>
    <template v-if="display==='visual'&&rows.length"><GeneVisual v-if="visualIndex" :index="visualIndex" :colors="data.kinetics_colors" :level="mode" :dataset="dataset" :gene="gene" :branch="branch" :parameter="parameter" :rows="rows"/><p v-else-if="visualLoading" role="status">Loading visual…</p><p v-else role="alert">{{ visualError||'No visual is available for this gene.' }}</p></template>
    <p v-if="!rows.length" class="empty">No methods match these selections.</p>
    <details v-if="mode==='kinetics'&&scope==='expanded'" class="source-details"><summary>Gene coverage</summary><p class="small">Counts below describe each method's available genes, including when scores above use a fixed common-gene set.</p><div class="table-scroll"><table><thead><tr><th>Method</th><th>Candidate</th><th>Present</th><th>Valid</th><th>Undefined</th><th>Unavailable</th></tr></thead><tbody><tr v-for="r in coverageRows" :key="r.method"><td>{{ names[r.method] }}</td><td>{{ r.candidate_n }}</td><td>{{ r.present_n }}</td><td>{{ r.n_genes }}</td><td>{{ r.undefined_n }}</td><td>{{ r.unavailable_n }}</td></tr></tbody></table></div></details>
    <p class="reference-link"><a href="#documentation">Metric definitions & evaluation settings</a> · <a href="#downloads">Download result tables</a></p>
  </section></div>
  <dialog v-if="enlarged" ref="dialog" class="image-dialog" :aria-label="`${enlarged.name} streamlines`" @cancel.prevent="enlarged=null" @click.self="enlarged=null"><div><button autofocus class="button secondary" @click="enlarged=null">Close</button><h3>{{ datasetName }} · {{ enlarged.name }}</h3><img :src="`./${enlarged.asset.path}`" :alt="`${enlarged.name} streamlines`"></div></dialog>
</template>
