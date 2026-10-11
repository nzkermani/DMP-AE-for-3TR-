---
name: interpret-pathway-enrichment
description: Pathway-level analysis for longitudinal multi-omics — per-sample ssGSEA/GSVA scores saved as a new analysable layer (then run through limma and LongOmicsTraj like any feature layer), plus preranked fgsea on every longitudinal contrast. Works for genes and proteins with MSigDB or any GMT, and for metabolites with a user-supplied GMT. Use when the user asks which pathways or biological processes change over visits or differ between arms.
---

# Pathway enrichment for repeated-visit data

```bash
Rscript .claude/skills/interpret-pathway-enrichment/scripts/pathway_scores.R \
    study_config.json <layer> <gene_sets> [ssgsea|gsva] [min_size=10] [max_size=500]
```
`gene_sets` = path to a `.gmt` file, or `msigdbr:<collection>` (e.g.
`msigdbr:H` hallmarks, `msigdbr:C2`). The script picks the right `msigdbr`
argument name for old (`category`) and new (`collection`) versions and the
right GSVA API (`ssgseaParam`/`gsvaParam` for GSVA ≥ 1.50, else the legacy call).

Outputs: `results/data/pathways_<layer>_processed_se.rds`,
`results/pathways/<layer>/coverage.csv`, and `fgsea_<contrast>.csv` for every
DE contrast if `fgsea` is installed and the differential step was run.

Then analyse pathways longitudinally exactly as features:
```bash
Rscript .claude/skills/stats-longitudinal-de/scripts/longitudinal_de.R study_config.json pathways_<layer>
Rscript .claude/skills/cluster-trajectories-longomicstraj/scripts/run_longomicstraj.R study_config.json pathways_<layer> empirical,lmm
```

## Before running — ID matching
Feature IDs must match the gene-set members. Check `coverage.csv`: if the
median matched fraction is low (< 0.5 for transcriptomics), the IDs are wrong
(Ensembl vs symbol, Olink/UniProt vs symbol). Ask the user for the mapping
table rather than guessing mappings.

Targeted proteomic panels (Olink 92/384, SomaScan subsets) cover a biased,
small fraction of each pathway: scores reflect the measured members only.
Say so when reporting, and keep `min_size` meaningful (≥ 5 matched).

Metabolites: no default database is bundled. Ask the user for a GMT built
from their annotation (KEGG / SMPDB / Reactome / HMDB classes) whose member
IDs match their metabolite IDs. With few annotated metabolites, consider
chemical-class sets instead of pathways.

## Method notes
- ssGSEA scores are z-scored per pathway across samples before saving, so
  `delta_min` in LongOmicsTraj is in SD units for pathway layers — use a
  larger value (e.g. 0.3) than for log2 features and state it.
- Single-sample scores are computed on all samples together; scores are not
  independent of the sample set — rerun after adding or removing samples.
- fgsea uses the limma moderated t as ranking statistic; seeded for
  reproducibility. Report NES and padj, and the leading edge for top hits.
