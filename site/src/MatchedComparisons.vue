<script setup>
import { ref } from 'vue'
import ImageViewer from './ImageViewer.vue'
const props = defineProps({ data: Object })
const results = props.data.comparisons
const enlarged = ref(null)
const files = './data/files/evaluation/'
const scores = files + 'cell/results/statistics/'
const gene = files + 'gene/results/crosslevel/'
const rates = files + 'kinetics_recovery/results/joint_recovery/'
const labels = {consistency:'Consistency Score', iccoh:'ICCoh', cbdir:'CBDir', transcosine:'TransCosine', batc:'BATC', gene_cbdir:'Gene CBDir', cell_cbdir:'Cell CBDir'}
const format = x => Number(x).toFixed(4)
const pvalue = x => Number(x) < .001 ? Number(x).toExponential(2) : Number(x).toPrecision(3)
</script>

<template>
  <div class="page-heading"><h1>Matched comparisons and statistics</h1></div>
  <section class="panel content-panel"><h2>Gene and cell directions · Figure S9</h2>
    <p>Across seven datasets, the nine configurations shared 35 scorable gene–dataset pairs. TIvelo led mean gene-level CBDir, while UniTVelo led mean cell-level CBDir and BATC. Both belonged to the three-metric Pareto front.</p>
    <button class="phase-image" aria-label="Enlarge Figure S9" @click="enlarged={...results.figure,title:'Figure S9 · Matched Gene and Cell outputs'}"><img :src="'./'+results.figure.path" alt="Matched Gene/Cell scores and Pareto-front resampling frequencies across seven real datasets"></button>
    <div class="guide-links"><a :href="gene+'dataset_scores.csv'" download>Matched scores</a><a :href="gene+'frontiers.csv'" download>Front frequencies</a><a :href="gene+'plot_data.csv'" download>Panel data</a><a :href="gene+'leave_one_out.csv'" download>Leave-one-out fronts</a><a v-if="data.downloads.figures.some(f=>f.path==='figures/suppFigNine.pdf')" href="./data/files/figures/suppFigNine.pdf" download>Figure S9 PDF</a></div>
  </section>
  <section class="panel content-panel"><h2>Overall differences · Table S10</h2>
    <p>Datasets are paired blocks in each Friedman test. P-values use the Iman–Davenport F approximation and Holm adjustment within each panel. W measures agreement in ranks.</p>
    <div v-for="(tests,i) in results.omnibus" :key="i"><h3>{{i===0?'A · Cell scores, 14 datasets and 20 configurations':'B · Matched Gene/Cell scores, seven datasets and nine configurations'}}</h3>
      <div class="table-scroll"><table><thead><tr><th>Metric</th><th class="numeric">Friedman Q</th><th class="numeric">Kendall W</th><th class="numeric">Holm-adjusted p</th></tr></thead><tbody><tr v-for="r in tests" :key="r.metric"><td>{{labels[r.metric]}}</td><td class="numeric">{{format(r.friedman_q)}}</td><td class="numeric">{{format(r.kendall_w)}}</td><td class="numeric">{{pvalue(r.p_holm)}}</td></tr></tbody></table></div>
    </div>
    <p>Across {{results.batc.comparisons}} leave-one-dataset-out comparisons, veloVI’s lowest rank by mean BATC was {{results.batc.velovi_worst}} and it ranked first {{results.batc.velovi_first}} times. Across individual datasets, {{results.batc.dataset_winners}} configurations ranked first at least once.</p>
    <div class="guide-links"><a :href="scores+'s10a.csv'" download>Table S10A</a><a :href="scores+'s10b.csv'" download>Table S10B</a><a :href="scores+'batc_dataset_ranks.csv'" download>Dataset ranks</a><a :href="scores+'batc_leave_one_out.csv'" download>Leave-one-out ranks</a></div>
  </section>
  <section class="panel content-panel"><h2>Jointly scorable transcription and degradation rates</h2>
    <p>Within each system or organoid branch, both methods have finite Spearman correlations for both rates on every included gene. Organoid branch means receive equal weight.</p>
    <div class="table-scroll"><table><thead><tr><th>System</th><th>Common genes</th><th>Method</th><th class="numeric">Mean α correlation</th><th class="numeric">Mean γ correlation</th></tr></thead><tbody>
      <template v-for="(system,key) in results.kinetics" :key="key"><tr v-for="(means,method) in system.means" :key="method"><td>{{key==='rpe1'?'RPE1':'Organoid'}}</td><td>{{system.components.length===1?system.components[0].common_genes:system.components.map(c=>'branch '+c.branch+': '+c.common_genes).join(', ')}}</td><td>{{method}}</td><td class="numeric">{{format(means.alpha)}}</td><td class="numeric">{{format(means.gamma)}}</td></tr></template>
    </tbody></table></div><div class="guide-links"><a :href="rates+'common_gene_scores.csv'" download>Common-gene scores</a><a :href="rates+'common_gene_comparison.json'" download>System summaries</a></div>
  </section>
  <ImageViewer v-if="enlarged" :image="enlarged" @close="enlarged=null"/>
</template>
