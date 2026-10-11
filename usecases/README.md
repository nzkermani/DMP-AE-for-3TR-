# Use cases — complete workflows

Each folder reproduces a published (or typical) analysis end to end with the
DMP-AE skills, on synthetic data that copy the study design and contain
known planted effects. A use case has:

- `study_config.json`, which points every skill at the use case's data and its own `results/`;
- `data/simulate_*.R`, a seeded generator, plus the generated CSVs and `planted_truth.csv`;
- an R Markdown notebook that calls the skill scripts in the paper's order,
  shows the outputs next to the paper's tables and figures, and ends with a
  replication scorecard;
- `README.md`, which maps each paper step to a skill and lists what had to be added.

| # | folder | paper | new skills it added |
|---|---|---|---|
| 01 | `01_breathomics_biologics` | Vassileva *et al.*, *Respir Med* 2026 — exhaled VOCs before/after benralizumab and dupilumab | `prep-breathomics-gcms`, `ml-multilevel-splsda`, `stats-clinical-association` |

## Adding a use case
1. Make `usecases/<nn>_<name>/` with a config whose `output_dir` points inside the folder.
2. Write a seeded simulator that matches the paper's design: n, visits,
   arms, layers, batch, missing visits. Plant the reported effects and
   write `planted_truth.csv`.
3. Map each method in the paper to a skill. When a method is missing, add a
   skill under `.claude/skills/` (SKILL.md + script). Put shared code in `R/`
   and an MCP tool in `mcp/longomics_tools.R` if Claude needs to query its
   results, then register the skill in `omics-analysis-guide`, `CLAUDE.md` and `README.md`.
4. Write the notebook with the `skill()` helper from use case 01, so cached
   steps are skipped. End it with a scorecard of paper vs this run.
5. Add a fast check to `tests/run_all_tests.R`.
