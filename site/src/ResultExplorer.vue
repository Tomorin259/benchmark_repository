<script setup>
import { computed, nextTick, onMounted, ref, watch } from 'vue'
import GeneVisual from './GeneVisual.vue'
import OverallResults from './OverallResults.vue'
import CrossLevel from './CrossLevel.vue'
import StrategyResults from './StrategyResults.vue'
import KineticsStreams from './KineticsStreams.vue'

const props = defineProps({ data:Object, initialLevel:{type:String,default:'cell'} })
const levels=['overall','cell','gene','kinetics','cross']
const mode=ref(levels.includes(props.initialLevel)?props.initialLevel:'cell')
const dataset=ref(mode.value==='kinetics'?'sceu_rpe1':'pancreas'), metric=ref(mode.value==='kinetics'?'score':'cbdir')
const scope=ref('full'), comparison=ref('available'), branch=ref(''), parameter=ref('alpha'), gene=ref('__summary__')
const query=ref(''), order=ref('method'), display=ref('table'), selected=ref(props.data.methods.map(m=>m.config_id).concat('dynamo_labelled'))
const visualIndex=ref(null), visualLoading=ref(false), visualError=ref(''), ready=ref(false)
const streamMethods=ref([]), compareOpen=ref(false), enlarged=ref(null), dialog=ref(null), copied=ref(false)
const summaryTask=ref('cell'), summaryStatistic=ref('mean'), summaryView=ref('table'), crossCase=ref(props.data.cross_cases[0].id)
const allReal=computed(()=>mode.value==='cell'&&dataset.value==='__all_real__')
const standard=computed(()=>['cell','gene','kinetics'].includes(mode.value))
const specialKinetics=computed(()=>mode.value==='kinetics'&&['strategy','streams'].includes(scope.value))
const strategyQuery=ref(''), strategyClass=ref(''), kineticMethods=ref([]), kineticCompareOpen=ref(false)
const names=Object.fromEntries(props.data.methods.map(m=>[m.config_id,m.manuscript_name]))
names.dynamo_labelled='Dynamo (scEU-seq)'
const kineticChoices=computed(()=>props.data.kinetics_streams.filter(r=>r.dataset===dataset.value).sort((a,b)=>props.data.kinetics_method_order.indexOf(a.method)-props.data.kinetics_method_order.indexOf(b.method)))
const kineticPreview=computed({get:()=>kineticMethods.value[0]??'',set:value=>{kineticMethods.value=[value];kineticCompareOpen.value=false}})
const localCase=computed(()=>mode.value==='gene'?props.data.real_cases.find(c=>c.dataset===dataset.value&&c.gene===gene.value):undefined)
watch(kineticChoices,choices=>{const valid=kineticMethods.value.filter(m=>choices.some(r=>r.method===m));kineticMethods.value=valid.length?valid:choices.slice(0,1).map(r=>r.method)},{immediate:true})
watch([mode,dataset,scope],()=>{if(mode.value==='kinetics'&&scope.value==='strategy'&&dataset.value!=='sceu_rpe1')scope.value='full'})
const cellMetrics=[['consistency','Consistency'],['iccoh','ICCoh'],['cbdir','CBDir'],['transcosine','TransCosine'],['batc','BATC']]
const comparisonNames={available:'Each method’s available genes',core_conventional:'Fixed common set · core RNA methods',all_conventional:'Fixed common set · all RNA methods',all_displayed_methods:'Fixed common set · RNA methods + scEU comparator'}
const activeMetric=computed(()=>props.data.content.metrics[mode.value].find(m=>m.key===metric.value))
const tables=computed(()=>mode.value==='cell'?[...props.data.cell,...props.data.cell_simulation]:(props.data[mode.value]??[]))
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
  query.value=''; gene.value='__summary__'; display.value='table'; streamMethods.value=[]; compareOpen.value=false
})
watch(()=>props.initialLevel,value=>{if(levels.includes(value)&&value!==mode.value)mode.value=value})
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
const scopeLabel=computed(()=>mode.value==='cell'?'':mode.value==='gene'?metric.value==='cossim'?'Simulation ground truth':gene.value==='__summary__'?'Mean · Finite evaluated genes':'Single gene · Reference 30-NN':scope.value==='selected'?`Selected examples · ${dataset.value==='sceu_rpe1'?19:branch.value==='1'?6:8} genes`:scope.value==='full'?'Canonical scores':`Extended reference · ${props.data.kinetics_expanded_scope[dataset.value]} candidate genes`)
const diagnosticFor=id=>Object.keys(props.data.diagnostics).find(key=>props.data.diagnostics[key].method===id)
const diagnosticLink=id=>`#diagnostics?method=${diagnosticFor(id)}&dataset=${dataset.value}`
const coverageRows=computed(()=>props.data.kinetics_summary.filter(r=>r.scope==='expanded'&&r.comparison==='available'&&r.dataset===dataset.value&&r.branch===branch.value&&r.parameter===parameter.value&&rows.value.some(x=>x.method===r.method)))
const streamRows=computed(()=>streamMethods.value.map(id=>rows.value.find(r=>r.method===id)).filter(Boolean).map(r=>({...r,asset:props.data.streams.find(s=>s.dataset===dataset.value&&s.method===r.method)})))
const previewChoices=computed(()=>props.data.methods.map(m=>unfilteredRows.value.find(r=>r.method===m.config_id)).filter(Boolean).map(r=>({...r,asset:props.data.streams.find(s=>s.dataset===dataset.value&&s.method===r.method)})).filter(r=>r.asset))
const previewMethod=computed({get:()=>streamMethods.value[0]??'',set:value=>{if(value){streamMethods.value=[value];compareOpen.value=false}}})
watch([mode,rows],([level,value])=>{
  const valid=streamMethods.value.filter(id=>value.some(r=>r.method===id))
  streamMethods.value=level==='cell'&&props.data.streams.some(s=>s.dataset===dataset.value)?(valid.length?valid:[value[0]?.method].filter(Boolean)):[]
},{immediate:true})
const format=value=>value===null||value===undefined?'—':value.toFixed(4)
const status=row=>metric.value==='cossim'?'':(row.status&&row.status!=='ok'?({gene_missing:'Gene absent from output',constant_prediction:'Constant predicted trajectory'}[row.status]||row.status.replaceAll('_',' ')):'')
const hasStatus=computed(()=>rows.value.some(r=>status(r)))
const bar=value=>({left:`${value<0?(value+1)*50:50}%`,width:`${Math.abs(value)*50}%`})
function chooseSubject(value){gene.value=value==='summary'?'__summary__':genes.value[0]??'';display.value=value==='summary'?'table':'visual'}
function sortMetric(key){if(metric.value===key)order.value=order.value==='descending'?'ascending':'descending';else{metric.value=key;order.value='descending'}}
async function preview(id){streamMethods.value=[id];compareOpen.value=false;display.value='table';await nextTick();document.getElementById('stream-preview')?.scrollIntoView({block:'start',behavior:'smooth'})}
function toggleCompare(id){if(streamMethods.value.includes(id)){if(streamMethods.value.length>1)streamMethods.value=streamMethods.value.filter(m=>m!==id);return}if(streamMethods.value.length<4)streamMethods.value=[...streamMethods.value,id]}
function toggleCompareMode(){if(compareOpen.value){streamMethods.value=streamMethods.value.slice(0,1);compareOpen.value=false}else compareOpen.value=true}
function toggleKineticCompare(){if(kineticCompareOpen.value){kineticMethods.value=kineticMethods.value.slice(0,1);kineticCompareOpen.value=false}else kineticCompareOpen.value=true}
function toggleKineticMethod(id){if(kineticMethods.value.includes(id)){if(kineticMethods.value.length>1)kineticMethods.value=kineticMethods.value.filter(m=>m!==id)}else if(kineticMethods.value.length<4)kineticMethods.value=[...kineticMethods.value,id]}
watch(enlarged,async value=>{await nextTick();if(value)dialog.value?.showModal()})
watch([display,gene],async()=>{
  if(display.value!=='visual'||gene.value==='__summary__'||visualIndex.value||visualLoading.value||!props.data.visuals)return
  visualLoading.value=true;visualError.value=''
  try{const response=await fetch(`./${props.data.visuals.index}`);if(!response.ok)throw Error(response.status);visualIndex.value=await response.json()}
  catch(e){console.error(e);visualError.value='The visual could not be loaded. Choose Table or try another gene.'}
  finally{visualLoading.value=false}
})
const shareUrl=computed(()=>{
  if(mode.value==='overall')return `${location.origin}${location.pathname}#results?${new URLSearchParams({level:'overall',task:summaryTask.value,statistic:summaryStatistic.value,view:summaryView.value})}`
  if(mode.value==='cross')return `${location.origin}${location.pathname}#results?${new URLSearchParams({level:'cross',case:crossCase.value})}`
  if(allReal.value)return `${location.origin}${location.pathname}#results?${new URLSearchParams({level:'cell',dataset:'__all_real__',statistic:summaryStatistic.value,view:summaryView.value})}`
  if(specialKinetics.value){
    const p=new URLSearchParams({level:'kinetics',dataset:dataset.value,scope:scope.value})
    if(scope.value==='strategy'){if(strategyQuery.value)p.set('geneSearch',strategyQuery.value);if(strategyClass.value)p.set('class',strategyClass.value)}
    else if(kineticMethods.value.length>1)p.set('compare',kineticMethods.value.join(','));else if(kineticMethods.value[0])p.set('method',kineticMethods.value[0])
    return `${location.origin}${location.pathname}#results?${p}`
  }
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
  if(['cell','gene','kinetics'].includes(p.get('task')))summaryTask.value=p.get('task')
  if(['mean','median'].includes(p.get('statistic')))summaryStatistic.value=p.get('statistic')
  if(['table','distribution'].includes(p.get('view')))summaryView.value=p.get('view')
  if([...props.data.cross_cases,...props.data.real_cases].some(c=>c.id===p.get('case')))crossCase.value=p.get('case')
  if(mode.value==='cell'&&p.get('dataset')==='__all_real__')dataset.value='__all_real__'
  if(datasets.value.some(d=>d.dataset_id===p.get('dataset')))dataset.value=p.get('dataset')
  if(['full','selected','expanded','strategy','streams'].includes(p.get('scope')))scope.value=p.get('scope')
  strategyQuery.value=p.get('geneSearch')??''
  if(props.data.kinetics_strategy.classes.includes(p.get('class')))strategyClass.value=p.get('class')
  await nextTick();if(branches.value.includes(p.get('branch')))branch.value=p.get('branch')
  await nextTick();if(parameters.value.includes(p.get('parameter')))parameter.value=p.get('parameter')
  if(comparisonNames[p.get('comparison')])comparison.value=p.get('comparison')
  if(metrics.value.some(m=>m[0]===p.get('metric')))metric.value=p.get('metric')
  await nextTick();if(genes.value.includes(p.get('gene'))||p.get('gene')==='__summary__'&&summaryAllowed.value)gene.value=p.get('gene')
  if(['method','ascending','descending'].includes(p.get('order')))order.value=p.get('order')
  if(mode.value!=='cell'&&p.has('methods'))selected.value=p.get('methods').split(',').filter(m=>names[m])
  query.value=mode.value==='cell'?'':p.get('search')??''
  await nextTick();if(['table','plot','visual'].includes(p.get('view'))&&!(gene.value==='__summary__'&&p.get('view')==='visual'))display.value=p.get('view')
  if(mode.value==='cell'){
    const requested=p.has('compare')?p.get('compare').split(','):p.has('method')?[p.get('method')]:p.has('streams')?p.get('streams').split(','):[]
    const valid=requested.filter(m=>rows.value.some(r=>r.method===m)).slice(0,4)
    if(valid.length){streamMethods.value=[...new Set(valid)];compareOpen.value=streamMethods.value.length>1}
  }
  if(mode.value==='kinetics'&&scope.value==='streams'){
    const requested=p.has('compare')?p.get('compare').split(','):p.has('method')?[p.get('method')]:[]
    const valid=[...new Set(requested.filter(m=>kineticChoices.value.some(r=>r.method===m)))].slice(0,4)
    if(valid.length){kineticMethods.value=valid;kineticCompareOpen.value=valid.length>1}
  }
  ready.value=true
  history.replaceState(null,'',shareUrl.value)
})
</script>

<template>
  <div class="page-heading"><h1>Results</h1><button class="button secondary" @click="copyLink">{{ copied?'Link copied':'Share this view' }}</button></div>
  <div class="tabs" aria-label="Evaluation tasks"><button v-for="[key,label] in [['cell','Cell-level'],['gene','Gene-level'],['kinetics','Kinetic-parameter recovery'],['cross','Cross-level'],['overall','Overall']]" :key="key" :class="{active:mode===key}" :aria-pressed="mode===key" @click="mode=key">{{ label }}</button></div>
  <OverallResults v-if="mode==='overall'" :data="data" v-model:task="summaryTask" v-model:statistic="summaryStatistic" v-model:view="summaryView"/>
  <CrossLevel v-else-if="mode==='cross'" :data="data" v-model:case-id="crossCase"/>
  <div v-else-if="standard" class="explorer"><aside class="filters panel">
    <label for="dataset">{{ mode==='kinetics'?'System':'Dataset' }}</label><select id="dataset" v-model="dataset"><option v-if="mode==='cell'" value="__all_real__">All real datasets</option><template v-if="mode==='kinetics'"><option v-for="d in datasets" :key="d.dataset_id" :value="d.dataset_id">{{d.paper_name}}</option></template><template v-else><optgroup v-for="[label,simulated] in [['Real datasets',false],['Simulations',true]]" :key="label" :label="label"><option v-for="d in datasets.filter(d=>d.dataset_id.startsWith('sim_')===simulated)" :key="d.dataset_id" :value="d.dataset_id">{{d.paper_name}}</option></optgroup></template></select>
    <template v-if="mode==='kinetics'"><label for="scope">Analysis</label><select id="scope" v-model="scope"><option value="full">Canonical benchmark scores</option><option value="selected">Selected example genes</option><option value="expanded">Extended reference-gene analysis</option><option v-if="dataset==='sceu_rpe1'" value="strategy">Turnover strategy</option><option value="streams">Velocity fields</option></select>
      <template v-if="!specialKinetics">
      <template v-if="scope==='expanded'"><label for="comparison">Gene comparison</label><select id="comparison" v-model="comparison"><option v-for="(label,key) in comparisonNames" :key="key" :value="key">{{ label }}</option></select></template>
      <template v-if="dataset==='sceu_organoid'"><label for="branch">Differentiation branch</label><select id="branch" v-model="branch"><option value="1">Branch 1 · Secretory</option><option value="2">Branch 2 · Enterocyte</option></select></template>
      <label for="parameter">Parameter</label><select id="parameter" v-model="parameter"><option v-for="p in parameters" :key="p" :value="p">{{ p==='alpha'?'α · Transcription':'γ · Degradation' }}</option></select>
      </template>
    </template>
    <template v-if="mode!=='cell'&&!specialKinetics"><label>Evaluation object</label><div class="subject-switch"><button :disabled="!summaryAllowed" :aria-pressed="gene==='__summary__'" @click="chooseSubject('summary')">Summary</button><button :aria-pressed="gene!=='__summary__'" @click="chooseSubject('single')">Single gene</button></div><template v-if="gene!=='__summary__'"><label for="gene">Gene</label><select id="gene" v-model="gene"><option v-for="g in genes" :key="g">{{ g }}</option></select></template></template>
    <template v-if="mode==='kinetics'&&scope==='streams'"><label for="kinetic-preview-method">Preview method</label><select id="kinetic-preview-method" v-model="kineticPreview"><option v-for="r in kineticChoices" :key="r.method" :value="r.method">{{names[r.method]}}</option></select><button class="text-button compare-toggle" :aria-expanded="kineticCompareOpen" aria-controls="kinetic-compare-picker" @click="toggleKineticCompare">{{kineticCompareOpen?'Close comparison':'Compare methods'}}</button><fieldset v-if="kineticCompareOpen" id="kinetic-compare-picker" class="compare-picker side-compare-picker"><legend>Methods to show (up to four)</legend><label v-for="r in kineticChoices" :key="r.method"><input type="checkbox" :checked="kineticMethods.includes(r.method)" :disabled="kineticMethods.length>=4&&!kineticMethods.includes(r.method)" @change="toggleKineticMethod(r.method)">{{names[r.method]}}</label></fieldset></template>
    <template v-if="mode==='cell'&&previewChoices.length"><label for="preview-method">Preview method</label><select id="preview-method" v-model="previewMethod"><option v-for="r in previewChoices" :key="r.method" :value="r.method">{{ r.name }}</option></select><button class="text-button compare-toggle" :aria-expanded="compareOpen" aria-controls="cell-compare-picker" @click="toggleCompareMode">{{ compareOpen?'Close comparison':'Compare methods' }}</button><fieldset v-if="compareOpen" id="cell-compare-picker" class="compare-picker side-compare-picker"><legend>Methods to show (up to four)</legend><label v-for="r in previewChoices" :key="r.method"><input type="checkbox" :checked="streamMethods.includes(r.method)" :disabled="streamMethods.length>=4&&!streamMethods.includes(r.method)" @change="toggleCompare(r.method)">{{ r.name }}</label></fieldset></template>
    <template v-if="mode==='gene'"><label for="metric">Metric</label><select id="metric" v-model="metric"><option v-for="m in metrics" :key="m[0]" :value="m[0]">{{ m[1] }}</option></select></template>
    <template v-if="mode!=='cell'&&!specialKinetics"><label for="method-search">Search methods</label><input id="method-search" v-model="query" placeholder="Method name">
      <details class="method-select"><summary>Methods <span>{{ rows.length }} / {{ availableMethods.length }}</span></summary><div class="method-actions"><button @click="selected=[...availableMethods]">Select all</button><button @click="selected=[]">Clear</button></div><label v-for="m in availableMethods" :key="m" class="checkbox"><input v-model="selected" type="checkbox" :value="m">{{ names[m] }}</label></details>
    </template>
  </aside>
  <StrategyResults v-if="mode==='kinetics'&&scope==='strategy'" :data="data" v-model:query="strategyQuery" v-model:reference-class="strategyClass"/>
  <KineticsStreams v-else-if="mode==='kinetics'&&scope==='streams'" :data="data" :dataset="dataset" :methods="kineticMethods"/>
  <OverallResults v-else-if="allReal" :data="data" cell-only v-model:statistic="summaryStatistic" v-model:view="summaryView"/>
  <section v-else class="panel result-panel"><div class="result-heading"><div><h2>{{ datasetName }}<span v-if="mode!=='cell'&&gene!=='__summary__'"> · {{ gene }}</span></h2><div class="scope-tags"><span v-if="scopeLabel" class="tag">{{ scopeLabel }}</span><span v-if="mode==='cell'&&dataset.startsWith('sim_')" class="tag">Simulation · Matched method outputs</span><span v-if="mode==='kinetics'&&branch" class="tag">{{ branch==='1'?'Secretory':'Enterocyte' }} branch</span></div></div><div class="guide-links"><a v-if="localCase" class="text-link" :href="`#results?level=cross&case=${localCase.id}`">Local comparison →</a><a class="text-link" href="#directory">Data & methods →</a></div></div>
    <div v-if="mode!=='cell'" class="result-context"><details class="metric-help" @keydown.esc="$event.currentTarget.open=false"><summary>{{ activeMetric.label }} ⓘ</summary><p>{{ activeMetric.definition }}</p><a href="#documentation">Metric reference →</a></details><span v-if="mode==='kinetics'&&scope==='expanded'&&comparison!=='available'" class="tag" :title="'Fixed for ' + commonMethods.map(m=>names[m]).join(', ')">Fixed common genes · {{ commonGenes.length }}</span></div>
    <section v-if="mode==='cell'&&streamRows.length" id="stream-preview" class="stream-preview">
      <p v-if="dataset==='forebrain'" class="small">Display: PC1/PC2 · CBDir: 30 PCs · BATC: 2D UMAP</p>
      <div class="stream-grid" :class="{single:streamRows.length===1}"><article v-for="r in streamRows" :key="r.method"><div class="image-heading"><h3>{{ r.name }}</h3><a v-if="diagnosticFor(r.method)" class="small" :href="diagnosticLink(r.method)">Diagnostics →</a></div><button v-if="r.asset" class="phase-image" :aria-label="`Enlarge ${r.name} streamlines`" @click="enlarged=r"><img :src="`./${r.asset.path}`" :alt="`${datasetName}: ${r.name} cell velocity field`" loading="lazy"></button><p v-else>Streamlines unavailable.</p><dl class="metric-values"><div v-for="[key,label] in cellMetrics" :key="key"><dt>{{ label }}</dt><dd>{{ format(r[key]) }}</dd></div></dl><details v-if="r.n_cells<r.n_reference" class="coverage-note"><summary>Partial cell coverage</summary><p>{{ r.n_cells }} of {{ r.n_reference }} reference cells are aligned for this method. This is cell coverage, not the valid-sample count of each metric.</p></details></article></div>
      <div class="cell-legend" v-if="streamRows[0]?.asset"><span v-for="item in streamRows[0].asset.legend" :key="item.label"><i :style="{background:item.color}"></i>{{ item.label }}</span></div>
    </section>
    <div class="table-toolbar"><div class="view-switch"><button class="text-button" :aria-pressed="display==='table'" @click="display='table'">{{ mode==='cell'?'All metrics':'Table' }}</button><button class="text-button" :aria-pressed="display==='plot'" @click="display='plot'">Score plot</button><button v-if="mode!=='cell'" class="text-button" :disabled="gene==='__summary__'" :aria-pressed="display==='visual'" @click="display='visual'">{{ mode==='gene'?'Phase portraits':'Trajectories' }}</button></div><div class="table-controls"><label v-if="mode==='cell'&&display==='plot'" class="order-label" for="plot-metric">Plot metric<select id="plot-metric" v-model="metric"><option v-for="[key,label] in cellMetrics" :key="key" :value="key">{{label}}</option></select></label><label class="order-label">Order<select v-model="order"><option value="method">Method</option><option value="descending">{{mode==='cell'?metricName+' ↓':'Score ↓'}}</option><option value="ascending">{{mode==='cell'?metricName+' ↑':'Score ↑'}}</option></select></label></div></div>
    <div v-if="display==='table'&&mode==='cell'" class="table-scroll"><table class="cell-scores"><thead><tr><th>Method</th><th class="numeric" v-for="[key,label] in cellMetrics" :key="key"><button @click="sortMetric(key)" :aria-label="`Sort by ${label}`">{{ label }}{{ metric===key&&order!=='method'?(order==='descending'?' ↓':' ↑'):'' }}</button></th></tr></thead><tbody><tr v-for="r in rows" :key="r.method"><td><button v-if="previewChoices.length" class="method-link" @click="preview(r.method)">{{ r.name }}</button><span v-else>{{r.name}}</span></td><td class="numeric" v-for="[key] in cellMetrics" :key="key">{{ format(r[key]) }}</td></tr></tbody></table></div>
    <div v-else-if="display==='table'" class="table-scroll"><table><thead><tr><th>Method</th><th class="numeric">{{ metricName }}</th><th v-if="metric!=='cossim'" class="numeric">{{ gene==='__summary__'?'Valid genes':mode==='gene'?'Source-cell observations':'Paired positions' }}</th><th v-if="hasStatus">Status</th></tr></thead><tbody><tr v-for="r in rows" :key="r.method"><td>{{ r.name }}<a v-if="diagnosticFor(r.method)" class="diagnostic-row-link" :href="diagnosticLink(r.method)">Diagnostics →</a><small v-if="r.method==='dynamo_labelled'" class="muted">Labeling-data comparator</small></td><td class="numeric" :title="String(r.value)">{{ format(r.value) }}</td><td v-if="metric!=='cossim'" class="numeric">{{ r.n_genes??r.n_cells??r.n_samples??'—' }}</td><td v-if="hasStatus">{{ status(r) }}</td></tr></tbody></table></div>
    <div v-if="display==='plot'&&rows.length" class="comparison"><div class="axis"><span>−1</span><span>0</span><span>1</span></div><div class="chart-row" v-for="r in rows" :key="r.method"><span class="chart-label">{{ r.name }}</span><div class="track"><div class="zero"></div><div v-if="r.value!==null" class="score-bar" :class="{negative:r.value<0}" :style="bar(r.value)"></div></div><span class="score-number">{{ format(r.value) }}</span></div></div>
    <template v-if="display==='visual'&&rows.length"><GeneVisual v-if="visualIndex" :index="visualIndex" :colors="data.kinetics_colors" :level="mode" :dataset="dataset" :gene="gene" :branch="branch" :parameter="parameter" :rows="rows"/><p v-else-if="visualLoading" role="status">Loading visual…</p><p v-else role="alert">{{ visualError||'No visual is available for this gene.' }}</p></template>
    <p v-if="!rows.length" class="empty">No methods match these selections.</p>
    <details v-if="mode==='kinetics'&&scope==='expanded'" class="source-details"><summary>Gene coverage</summary><span class="small muted">Available genes per method</span><div class="table-scroll"><table><thead><tr><th>Method</th><th>Candidate</th><th>Present</th><th>Valid</th><th>Undefined</th><th>Unavailable</th></tr></thead><tbody><tr v-for="r in coverageRows" :key="r.method"><td>{{ names[r.method] }}</td><td>{{ r.candidate_n }}</td><td>{{ r.present_n }}</td><td>{{ r.n_genes }}</td><td>{{ r.undefined_n }}</td><td>{{ r.unavailable_n }}</td></tr></tbody></table></div></details>
    <p class="reference-link"><a href="#documentation">Metric definitions & evaluation settings</a> · <a href="#downloads">Download result tables</a></p>
  </section></div>
  <dialog v-if="enlarged" ref="dialog" class="image-dialog" :aria-label="`${enlarged.name} streamlines`" @cancel.prevent="enlarged=null" @click.self="enlarged=null"><div><button autofocus class="button secondary" @click="enlarged=null">Close</button><h3>{{ datasetName }} · {{ enlarged.name }}</h3><img :src="`./${enlarged.asset.path}`" :alt="`${enlarged.name} streamlines`"></div></dialog>
</template>
