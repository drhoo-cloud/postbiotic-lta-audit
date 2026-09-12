# postbiotic-lta-audit

An audit of what the public sequence record can and cannot say about
lipoteichoic acid D-alanylation, and the analysis behind

> Nguyen TTM, Kim H, Kim H, Chung DK, Yi T-H.
> *Standardizing postbiotics: a genome cannot supply the value,
> and lipoteichoic acid shows what to state instead.* Submitted to Trends in Food Science & Technology.

Every number in the manuscript and its supplement is reproducible from this
repository. No new sequence and no new laboratory measurement were generated:
the analysis runs on RefSeq assemblies that were already deposited and on
substitution values that were already published.

## Layout

    config/               species list and every calling threshold
    scripts/              the screen, in the order it runs
    scripts/literature/   the PubMed sweep and lipoprotein predictor
    figures/              one script per panel
    data/tables/          the supplementary tables the manuscript cites

## What the screen asks

Not what the D-alanine substitution level of a strain is — that is a
physiological state and the sequence does not carry it — but whether a strain
has any means of removing D-alanine at all. The question divides into three
tests: is the removal enzyme present, can it reach its substrate, and how
conserved are the residues that hold that substrate.

## Running it

    pip install -r requirements.txt
    python scripts/01_collect_accessions.py
    python scripts/02_layer1_presence.py       # resumable; writes data/layer1.jsonl
    python scripts/03_layer2_localization.py   # set runtime.contact_email first
    python scripts/04_layer3_contacts.py
    python scripts/05_aggregate.py             # writes Tables S9a-b
    python scripts/06_dispersion_test.py       # writes Table S11d

Steps 1 to 5 take several hours and need network access to NCBI and EBI.
Step 6 does not: it runs in seconds from `data/tables/TableS1_structure_map.csv`
and `figures/fig2_data.json`, both included here, so the test reported in
Section 3 of the manuscript can be checked without rerunning anything.

## Figures

    python figures/make_figures.py            # everything, in order

or one at a time:

    python figures/make_fig2_data.py          # regenerates fig2_data.json from Table S9d
    python figures/fig1a_ester_chemistry.py   # Figure 1A
    python figures/fig1b_five_strains.py      # Figure 1B
    python figures/fig2_localization_and_process.py   # Figure 2
    python figures/fig3_measured_vs_genomic.py        # Figure 3
    python figures/figS1_score_distribution.py        # Figure S1 (needs Table S9a)
    python figures/figS2_contact_residues.py          # Figure S2
    python figures/figS3_mprf.py                      # Figure S3

Each writes a 600 dpi PNG and an LZW-compressed TIFF at the width used in the
manuscript. Figure 1 is composed from its two panel scripts by make_figures.py.

## Calling criteria and their limits

Every cut-off lives in `config/thresholds.yaml`. Nothing is hard-coded in the
scripts; if you find a hard-coded threshold, it is a bug, and the commit history
of that file is the record of what was called and when.

Two limits are stated rather than buried. The criteria were derived from a pilot
scan of 110 complete genomes, so they were fixed before the expanded scan but not
before any genome had been seen; the pre-registration record, with the eight prior
predictions and what became of each, is Methods S2 of the manuscript supplement.
And the localization call is a prediction from sequence topology, not a
measurement: Section 6 of the manuscript states the experiment that would test it.

## Data provenance

`data/tables/TableS1_structure_map.csv` compiles published measurements of
lipoteichoic acid structure; each row carries its primary source. The other
tables are outputs of the screen. `figures/fig2_data.json` is derived from
`data/tables/TableS9d.csv` and is regenerated, never edited by hand.

## Citing this repository

The release that accompanies the submitted manuscript is archived on Zenodo.
Cite the DOI rather than this URL, so that the version you used is the version
a reader gets.

    https://doi.org/10.5281/zenodo.22719721

That DOI resolves to the most recent version. The v1.0.0 release that
accompanies the submitted manuscript is 10.5281/zenodo.22719722.

`CITATION.cff` above carries the same information in machine-readable form;
GitHub renders it under "Cite this repository".

## Licence

MIT, see `LICENSE`. The supplementary tables under `data/tables/` are
compilations from the published literature; the sources for each row are given
in the manuscript's supplementary material.
