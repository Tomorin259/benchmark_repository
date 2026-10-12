<script setup>
import { computed, nextTick, ref, watch } from 'vue'

const props = defineProps({ index: Object, level: String, dataset: String, gene: String,
  branch: String, parameter: String, rows: Array, colors: Object })
const curve = ref(null)
const error = ref('')
const loading = ref(false)
const from = ref(0)
const to = ref(1)
const pointer = ref(null)
const dragStart = ref(null)
const enlarged = ref(null)
const dialog = ref(null)
watch(enlarged, async value => { await nextTick(); if (value) dialog.value?.showModal() })
const format = value => value === null || value === undefined ? '—' : value.toFixed(4)
const phases = computed(() => props.rows.map(row => ({ ...row,
  asset: props.index.phases.find(p => p.dataset === props.dataset && p.gene === props.gene && p.method === row.method) })))
const binding = computed(() => props.index.curves.find(c => c.dataset === props.dataset && c.gene === props.gene &&
  c.branch === props.branch && c.parameter === props.parameter))
watch(binding, async (value, previous, cleanup) => {
  const controller = new AbortController()
  cleanup(() => controller.abort())
  curve.value = null
  error.value = ''
  from.value = 0
  to.value = 1
  pointer.value = null
  loading.value = !!value
  if (!value) return
  try {
    const response = await fetch(`./${value.path}`, { signal: controller.signal })
    if (!response.ok) throw new Error(`HTTP ${response.status}`)
    curve.value = await response.json()
  } catch (e) {
    if (e.name !== 'AbortError') error.value = e.message
  } finally { if (!controller.signal.aborted) loading.value = false }
}, { immediate: true })
const series = computed(() => {
  const saved = curve.value?.series ?? []
  return [saved.find(s => s.method === 'reference'),
    ...props.rows.map(row => saved.find(s => s.method === row.method))].filter(Boolean).map(s => ({ ...s,
    name: s.method === 'reference' ? 'scEU-seq reference' : props.rows.find(r => r.method === s.method)?.name,
    color: s.method === 'reference' ? '#111827' : props.colors[s.method],
    score: props.rows.find(r => r.method === s.method)?.value }))
})
const xRange = computed(() => {
  const lower = Math.min(0.99, Math.max(0, Number(from.value) || 0))
  return [lower, Math.min(1, Math.max(lower + 0.01, Number(to.value) || 1))]
})
const separateAxes = computed(() => props.dataset === 'sceu_organoid')
const plotRight = computed(() => separateAxes.value ? 720 : 770)
function bounds(items) {
  const values = items.flatMap(s => s.y.filter((v, i) => v !== null &&
    s.x[i] >= xRange.value[0] && s.x[i] <= xRange.value[1]))
  if (!values.length) return [-1, 1]
  const low = Math.min(...values), high = Math.max(...values)
  const padding = (high - low || Math.abs(low) || 1) * 0.08
  return [low - padding, high + padding]
}
const yRange = computed(() => bounds(separateAxes.value ? series.value.filter(s => s.method !== 'reference') : series.value))
const referenceRange = computed(() => bounds(series.value.filter(s => s.method === 'reference')))
const x = value => 75 + (value - xRange.value[0]) / (xRange.value[1] - xRange.value[0]) * (plotRight.value - 75)
const y = (value, reference = false) => {
  const range = reference && separateAxes.value ? referenceRange.value : yRange.value
  return 360 - (value - range[0]) / (range[1] - range[0]) * 305
}
const xTicks = computed(() => Array.from({ length: 6 }, (_, i) => xRange.value[0] + i / 5 * (xRange.value[1] - xRange.value[0])))
const yTicks = computed(() => Array.from({ length: 5 }, (_, i) => yRange.value[0] + i / 4 * (yRange.value[1] - yRange.value[0])))
const referenceTicks = computed(() => Array.from({ length: 5 }, (_, i) => referenceRange.value[0] + i / 4 * (referenceRange.value[1] - referenceRange.value[0])))
function path(series) {
  let active = false
  return series.y.map((v, i) => {
    if (v === null || series.x[i] < xRange.value[0] || series.x[i] > xRange.value[1]) { active = false; return '' }
    const result = `${active ? 'L' : 'M'}${x(series.x[i]).toFixed(3)},${y(v, series.method === 'reference').toFixed(3)}`
    active = true
    return result
  }).join(' ')
}
const hovered = computed(() => pointer.value === null ? [] : series.value.map(s => {
  const index = s.x.reduce((best, value, i) => Math.abs(value - pointer.value) < Math.abs(s.x[best] - pointer.value) ? i : best, 0)
  return { ...s, position: s.positions[index], relative: s.x[index], value: s.y[index] }
}))
function position(event) {
  const svg = event.currentTarget
  const point = svg.createSVGPoint()
  point.x = event.clientX
  point.y = event.clientY
  const local = point.matrixTransform(svg.getScreenCTM().inverse())
  return Math.max(xRange.value[0], Math.min(xRange.value[1], xRange.value[0] + (local.x - 75) / (plotRight.value - 75) * (xRange.value[1] - xRange.value[0])))
}
function finishZoom(event) {
  const end = position(event)
  if (dragStart.value !== null && Math.abs(end - dragStart.value) > 0.01) {
    from.value = Math.min(end, dragStart.value)
    to.value = Math.max(end, dragStart.value)
  }
  dragStart.value = null
}
</script>

<template>
  <div v-if="level === 'gene'" class="phase-grid">
    <article v-for="row in phases" :key="row.method" class="phase-card">
      <div class="phase-heading"><h3>{{ row.name }}</h3><span>{{ format(row.value) }}</span></div>
      <button v-if="row.asset?.path" class="phase-image" :aria-label="`Enlarge ${row.name} phase portrait`" @click="enlarged = row">
        <img :src="`./${row.asset.path}`" :alt="`${gene}: ${row.name} phase portrait`" loading="lazy" width="600" height="600">
      </button>
      <p v-else class="empty">{{ row.asset?.status === 'gene_absent' ? 'Gene absent in method output.' : 'Phase portrait unavailable.' }}</p>
    </article>
    <dialog v-if="enlarged" ref="dialog" class="image-dialog" :aria-label="`${enlarged.name} phase portrait`" @click.self="enlarged = null" @cancel.prevent="enlarged = null">
      <div><button autofocus class="button secondary" @click="enlarged = null">Close</button><h3>{{ gene }} · {{ enlarged.name }}</h3><img :src="`./${enlarged.asset.path}`" :alt="`${gene} phase portrait`"></div>
    </dialog>
  </div>
  <div v-else class="trajectory-view">
    <p v-if="loading" role="status">Loading trajectories…</p>
    <p v-else-if="error" role="alert">Trajectory could not be loaded. Please select another gene or try again.</p>
    <p v-else-if="!curve">Trajectory unavailable for this gene.</p>
    <template v-else>
      <div v-if="!separateAxes" class="trajectory-controls"><span>{{ curve.transform }}</span></div>
      <div class="zoom-controls"><button class="text-button" @click="from = 0; to = 1">Reset zoom</button><span>Drag to zoom</span><details><summary>Set range</summary><label>From <input v-model.number="from" type="number" min="0" :max="to - 0.01" step="0.01" aria-label="Trajectory start" @change="from = xRange[0]; to = xRange[1]"></label><label>To <input v-model.number="to" type="number" :min="from + 0.01" max="1" step="0.01" aria-label="Trajectory end" @change="from = xRange[0]; to = xRange[1]"></label></details></div>
      <svg class="trajectory-chart" viewBox="0 0 810 455" role="img" :aria-label="`${gene} ${parameter} trajectories`" font-family="Arial, sans-serif" font-size="12" fill="#374151"
        @pointermove="pointer = position($event)" @pointerleave="pointer = null; dragStart = null" @pointerdown="dragStart = position($event)" @pointerup="finishZoom">
        <title>{{ dataset }} {{ branch ? `branch ${branch}` : '' }} · {{ gene }} · {{ parameter }} · {{ curve.transform }}</title>
        <text x="75" y="25" font-size="16" font-weight="600">{{ gene }} · {{ parameter === 'alpha' ? 'α' : 'γ' }}</text>
        <g v-for="tick in yTicks" :key="tick"><line x1="75" :x2="plotRight" :y1="y(tick)" :y2="y(tick)" stroke="#e5e7eb"/><text x="65" :y="y(tick) + 4" text-anchor="end">{{ tick.toFixed(2) }}</text></g>
        <g v-if="separateAxes"><line :x1="plotRight" :x2="plotRight" y1="55" y2="360" stroke="#111827"/><text v-for="tick in referenceTicks" :key="tick" :x="plotRight + 10" :y="y(tick,true) + 4">{{ tick.toFixed(2) }}</text><text transform="translate(799,208) rotate(-90)" text-anchor="middle">Reference (Table S4)</text></g>
        <g v-for="tick in xTicks" :key="tick"><text :x="x(tick)" y="380" text-anchor="middle">{{ tick.toFixed(2) }}</text></g>
        <line x1="75" :x2="plotRight" y1="360" y2="360" stroke="#9ca3af"/>
        <text x="423" y="401" text-anchor="middle">Relative position</text>
        <text transform="translate(18,208) rotate(-90)" text-anchor="middle">{{ dataset === 'sceu_rpe1' ? 'Display z-score' : 'Prediction (baseline ratio)' }}</text>
        <path v-for="s in series" :key="s.method" :d="path(s)" fill="none" :stroke="s.color" :stroke-width="s.method === 'reference' ? 2.8 : 1.8"/>
        <line v-if="pointer !== null" data-hover :x1="x(pointer)" :x2="x(pointer)" y1="55" y2="360" stroke="#6b7280" stroke-dasharray="4 4"/>
        <rect v-if="dragStart !== null && pointer !== null" data-hover :x="x(Math.min(dragStart,pointer))" y="55" :width="Math.abs(x(dragStart) - x(pointer))" height="305" fill="#2166ac" opacity="0.1"/>
        <g v-for="(s,i) in series" :key="`legend-${s.method}`" :transform="`translate(${75 + (i % 3) * 238},${423 + Math.floor(i / 3) * 21})`"><line x1="0" x2="17" stroke-width="2" :stroke="s.color"/><text x="23" y="4">{{ s.name }}</text></g>
      </svg>
      <div v-if="hovered.length" class="trajectory-tooltip"><span v-for="s in hovered" :key="s.method"><i :style="{ background: s.color }"></i>{{ s.name }}: {{ format(s.value) }} <small>(position {{ s.relative.toFixed(3) }})</small></span></div>
      <div class="trajectory-scores"><span v-for="s in series.filter(s => s.method !== 'reference')" :key="s.method"><i :style="{ background: s.color }"></i>{{ s.name }} · ρ {{ format(s.score) }}</span></div>
      <p v-if="rows.some(row => !series.some(s => s.method === row.method))" class="small muted">Some displayed methods have no saved trajectory for this gene.</p>
    </template>
  </div>
</template>
