<script setup>
import { ref, watch } from 'vue'
import ImageViewer from './ImageViewer.vue'

const props = defineProps({ reference: { type: Object, required: true }, task: String, open: Boolean })
const emit = defineEmits(['update:open'])
const enlarged = ref(false)
watch(() => [props.reference.dataset, props.task], () => { enlarged.value = false })
</script>

<template>
  <section class="transition-reference" aria-label="Annotated evaluation transitions">
    <div class="image-heading">
      <h3>Annotated transitions</h3>
      <button class="text-button" :aria-expanded="open" :aria-controls="'transition-body-' + reference.dataset" @click="emit('update:open', !open)">{{ open ? 'Hide' : 'Show' }}</button>
    </div>
    <div v-if="open" :id="'transition-body-' + reference.dataset">
      <button class="phase-image" :aria-label="`Enlarge ${reference.name} annotated transitions`" @click="enlarged = true">
        <img :src="'./' + reference.files.png.path" :alt="`${reference.name}: annotated evaluation transitions on the reference ${reference.frame.toUpperCase()}`" :width="reference.files.png.width" :height="reference.files.png.height">
      </button>
      <p class="transition-caption">Arrows indicate the annotated transitions used for {{ task === 'gene' ? 'CBDir-gene' : 'directional evaluation' }}.</p>
      <p v-if="reference.dataset.startsWith('sim_')" class="transition-caption">Node labels identify branch bins (A/B, 0–9), colored from early to late.</p>
      <div class="transition-tools">
        <details><summary>Transitions · {{ reference.edges.length }}</summary>
          <ul><li v-for="[source, target] in reference.edges" :key="source + '→' + target">{{ source }} → {{ target }}</li></ul>
        </details>
        <span class="transition-downloads"><a :href="'./' + reference.files.png.path" download>PNG</a><a :href="'./' + reference.files.svg.path" download>SVG</a></span>
      </div>
    </div>
    <ImageViewer v-if="enlarged" :image="{ title: reference.name + ' · Annotated transitions', path: reference.files.svg.path }" @close="enlarged = false"/>
  </section>
</template>
