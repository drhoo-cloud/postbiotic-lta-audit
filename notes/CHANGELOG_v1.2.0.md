# v1.2.0 — full re-run and supplementary renumbering (2026-09-23)

## What changed

- **Whole pipeline re-run** on 9,143 RefSeq assemblies (02 → 05, plus Phobius for all
  727 unique DltE proteins and the PDB 8AKH contact mapping). Every call in the previous
  release reproduced exactly: dltA, DltE at bit scores 90/110/130, long and short MprF,
  and all three localization classes. The only field that moved is `n_protein`, which
  changed in 279 genomes because RefSeq re-annotated them; no call depends on it.
- **Supplementary numbering changed** to run consecutively. Old `TableS3` → `TableS2`,
  old `TableS9a–e` → `TableS3a–e`, old `TableS11` → `TableS4`. File names follow.
- **`TableS3a_per_genome.csv` gains a `localization` column**, so the per-genome file now
  carries the layer-2 class alongside the layer-1 calls.
- **`TableS1_structure_map.csv` gains one row**: *L. delbrueckii* subsp. *lactis*
  ATCC 15808 as an n-butanol extract (69% D-alanine). The hot-phenol row for the same
  strain (21–27%) is unchanged. The pair is the only within-strain comparison of
  extraction methods in the compiled record.

## Numbers confirmed by the re-run

| quantity | value |
|---|---|
| lactic acid bacterial genomes | 7,404 (9,143 total, 1,739 controls) |
| DltE absent at bit score 110 | 2,464 (33.3%, Wilson 32.2–34.4) |
| DltE absent at bit score 90 / 130 | 2,006 (27.1%) / 3,608 (48.7%) |
| dlt operon absent | 107 (1.4%) |
| no extracytoplasmic removal | 4,188 (56.6%, Wilson 55.4–57.7) |
| unique DltE proteins | 727, carried by 5,857 genomes |
| controls returning DltE at 110 | 997 of 1,739 (positive control) |
| catalytic vs substrate-contact conservation | 99.9% vs 19.5% |

The previously published sensitivity figures of 28.2% and 49.8% are not reproduced by
any run and are superseded by 27.1% and 48.7% on the fixed denominator (all 7,404
lactic acid bacterial genomes, including the 107 without a dlt operon).
