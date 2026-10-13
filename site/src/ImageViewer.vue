<script setup>
import { onMounted, onUnmounted, ref } from 'vue'
defineProps({ image: { type: Object, required: true } })
const emit = defineEmits(['close'])
const dialog = ref(null)
const returnFocus = document.activeElement
onMounted(()=>dialog.value.showModal())
onUnmounted(()=>{if(returnFocus?.isConnected)returnFocus.focus()})
</script>

<template>
  <dialog ref="dialog" class="image-dialog" :aria-label="image.title" @cancel.prevent="emit('close')" @click.self="emit('close')"><div>
    <div class="image-toolbar"><h2>{{image.title}}</h2><a class="button secondary" :href="`./${image.path}`" download>Download</a><button autofocus class="button secondary" @click="emit('close')">Close</button></div>
    <img :src="`./${image.path}`" :alt="image.title">
  </div></dialog>
</template>
