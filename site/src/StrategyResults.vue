<script setup>
import { computed, ref, watch } from 'vue'
import ImageViewer from './ImageViewer.vue'

const props = defineProps({ data:Object, query:{type:String,default:''}, referenceClass:{type:String,default:''} })
const emit = defineEmits(['update:query','update:referenceClass'])
const strategy = props.data.kinetics_strategy
const names = Object.fromEntries(props.data.methods.map(m=>[m.config_id,m.manuscript_name]))
const enlarged = ref(null), page = ref(1)
const rows = computed(()=>strategy.genes.filter(r=>r.gene.toLowerCase().includes(props.query.toLowerCase())&&(!props.referenceClass||r.reference===props.referenceClass)))
const pageCount = computed(()=>Math.max(1,Math.ceil(rows.value.length/20)))
const displayed = computed(()=>rows.value.slice((page.value-1)*20,page.value*20))
watch([()=>props.query,()=>props.referenceClass],()=>{page.value=1})
const format = value=>value.toFixed(4)
</script>

<template>
  <section class="panel result-panel strategy-panel">
    <div class="result-heading"><h2>RPE1 · Turnover strategy</h2><span class="tag">{{strategy.n_genes}} matched genes</span></div>
    <div class="analysis-grid strategy-grid"><article v-for="m in strategy.methods" :key="m.method" class="strategy-method"><div class="image-heading"><h3>{{names[m.method]}}</h3><span class="small">Macro-F1 <strong>{{format(m.macro_f1)}}</strong></span></div><button class="phase-image" :aria-label="`Enlarge ${names[m.method]} strategy classification`" @click="enlarged={...m.image,title:names[m.method]+' · Turnover strategy'}"><img :src="'./'+m.image.path" :alt="`${names[m.method]}: Cooperative, Neutral and Destabilizing reference/prediction overlap`"></button><div class="table-scroll strategy-counts"><table><thead><tr><th>Class</th><th class="numeric">scEU only</th><th class="numeric"><abbr title="Prediction only">Pred. only</abbr></th><th class="numeric">Shared</th></tr></thead><tbody><tr v-for="c in m.counts" :key="c.class"><td>{{c.class}}</td><td class="numeric">{{c.reference_only}}</td><td class="numeric">{{c.prediction_only}}</td><td class="numeric">{{c.intersection}}</td></tr></tbody></table></div></article></div>
    <details class="source-details strategy-genes"><summary>Gene classification</summary>
      <div class="classification-filters"><label>Search genes<input :value="query" @input="emit('update:query',$event.target.value)" placeholder="Gene name"></label><label>Reference class<select :value="referenceClass" @change="emit('update:referenceClass',$event.target.value)"><option value="">All classes</option><option v-for="c in strategy.classes" :key="c">{{c}}</option></select></label></div>
      <div class="table-scroll"><table><thead><tr><th>Gene</th><th>scEU-seq reference</th><th v-for="m in strategy.methods" :key="m.method">{{names[m.method]}}</th></tr></thead><tbody><tr v-for="r in displayed" :key="r.gene"><td>{{r.gene}}</td><td>{{r.reference}}</td><td v-for="m in strategy.methods" :key="m.method" :class="{'class-mismatch':r.predictions[m.method]!==r.reference}">{{r.predictions[m.method]??'—'}}</td></tr></tbody></table></div>
      <p v-if="!rows.length" class="empty">No genes match these selections.</p>
      <div class="pagination"><span>{{rows.length}} {{rows.length===1?'gene':'genes'}}</span><button class="text-button" :disabled="page===1" @click="page--">Previous page</button><span>{{page}} / {{pageCount}}</span><button class="text-button" :disabled="page===pageCount" @click="page++">Next page</button></div>
    </details>
    <p class="reference-link"><a href="#downloads">Download classification tables →</a></p>
  </section>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>
