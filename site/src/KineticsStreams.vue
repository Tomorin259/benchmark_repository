<script setup>
import { computed, ref } from 'vue'
import ImageViewer from './ImageViewer.vue'
const props = defineProps({ data:Object, dataset:String, methods:Array })
const names = Object.fromEntries(props.data.methods.map(m=>[m.config_id,m.manuscript_name]))
names.dynamo_labelled='Dynamo (scEU-seq)'
const enlarged = ref(null)
const datasetName = computed(()=>props.data.datasets.find(d=>d.dataset_id===props.dataset)?.paper_name)
const rows = computed(()=>props.methods.map(method=>props.data.kinetics_streams.find(r=>r.dataset===props.dataset&&r.method===method)).filter(Boolean))
</script>

<template>
  <section class="panel result-panel kinetics-fields">
    <div class="result-heading"><h2>{{datasetName}}</h2><span class="tag">Velocity fields</span></div>
    <div class="stream-grid" :class="{single:rows.length===1}"><article v-for="r in rows" :key="r.method"><h3>{{names[r.method]}}</h3><button class="phase-image" :aria-label="`Enlarge ${names[r.method]} kinetic velocity field`" @click="enlarged={...r.image,title:datasetName+' · '+names[r.method]}"><img :src="'./'+r.image.path" :alt="`${datasetName}: ${names[r.method]} saved cell velocity field`" loading="lazy"></button></article></div>
  </section>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>
