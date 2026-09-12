#!/usr/bin/env python3
"""
queries.py -- the literature screen behind Table S11b.

Sweeps PubMed for reports of purified lipoteichoic acid in which at least two
strains are compared within one assay system. Its output is the entry pool
whose attrition Table S11b records.

Stage 2 of the screen. Stage 1 was a citation expansion from three anchor
papers; stage 2 is this method-vocabulary sweep, run against PubMed esearch.
Stage 3 re-expanded citations from every paper confirmed in stages 1-2.

Run:  python queries.py            # prints per-query totals
"""
import json, time, urllib.parse, urllib.request

ESEARCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
           "?db=pubmed&retmode=json&retmax=400&term=")

# ---- stage 1 : citation expansion anchors -------------------------------
ANCHORS = {
    "Morath 2001":  "11342591",   # J Exp Med 193:393
    "Morath 2002":  "12417633",   # J Exp Med 195:1635
    "Deininger 2003": "12646642", # J Immunol 170:4134
}


# ---- stage 2a : method-anchored sweep (14 queries) ----------------------
# Run first, to reach papers that do not cite the three anchors at all.
# This is the round that surfaced the 1988-2001 layer the citation network
# had missed entirely, including Bhakdi 1991.
QUERIES_METHOD = {
 # extraction and purification
 "butanol extraction":
   '("n-butanol" OR "butanol extraction") AND (lipoteichoic OR LTA)',
 "octyl-Sepharose":
   '("octyl-sepharose" OR "hydrophobic interaction chromatography") AND lipoteichoic',
 "phenol extraction":
   '("phenol extraction" OR "hot phenol") AND lipoteichoic',
 "HF / alkaline hydrolysis":
   '("hydrofluoric acid" OR "alkaline hydrolysis" OR "deacylated") AND lipoteichoic',
 # assay systems
 "HEK-Blue TLR2":
   '("HEK-Blue" OR "HEK293" OR "SEAP reporter") AND lipoteichoic',
 "CHO/CD14 TLR2":
   '"CHO/CD14" AND lipoteichoic',
 "whole blood":
   '"whole blood" AND lipoteichoic AND (cytokine OR IL-8 OR TNF)',
 "RAW264.7":
   '"RAW 264.7" AND lipoteichoic',
 "PBMC/monocyte":
   '(PBMC OR "mononuclear cells" OR monocyte) AND lipoteichoic AND (TNF OR IL-6)',
 # genetic controls
 "lgt mutant":
   '(lgt OR "lipoprotein diacylglyceryl transferase") AND lipoteichoic AND mutant',
 "ltaS mutant":
   '(ltaS OR "LTA synthase") AND mutant AND (cytokine OR TLR2 OR TNF)',
 "dlt mutant":
   '(dltA OR dltB OR dltC OR dltD OR "D-alanylation") AND mutant AND (TLR2 OR cytokine OR TNF)',
 # structural determination
 "NMR / MS structure":
   '(NMR OR "mass spectrometry" OR MALDI) AND lipoteichoic AND structure',
 # recent vocabulary
 "postbiotic":
   '(postbiotic OR paraprobiotic OR "cell-free supernatant") AND lipoteichoic',
}

# ---- stage 2b : readout, chemistry and genomics vocabulary (25 queries) --
QUERIES = {
 # readout systems -- anti-inflammatory and antagonism
 "anti-inflammatory":
   'lipoteichoic acid AND (anti-inflammatory OR "LPS-induced" OR tolerance OR "cross-tolerance")',
 "receptor antagonism":
   'lipoteichoic acid AND (antagonis* OR inhibit*) AND (TLR2 OR TLR4 OR "NF-kB")',
 "IL-10 regulatory":
   'lipoteichoic acid AND (IL-10 OR Treg OR regulatory OR immunomodulat*)',
 # readout systems -- mucosal immunity
 "mucosal/IgA":
   'lipoteichoic acid AND (IgA OR "Peyer" OR "lamina propria" OR mucosal OR "intestinal epithelial")',
 "dendritic cell":
   'lipoteichoic acid AND ("dendritic cell" OR "bone marrow-derived" OR BMDC OR "Langerhans")',
 "gut/colitis":
   'lipoteichoic acid AND (colitis OR DSS OR "gut barrier" OR "tight junction")',
 "skin/keratinocyte":
   'lipoteichoic acid AND (keratinocyte OR "atopic dermatitis" OR skin)',
 "osteoclast/bone":
   'lipoteichoic acid AND (osteoclast OR osteoblast OR "bone resorption")',
 "mast cell/allergy":
   'lipoteichoic acid AND ("mast cell" OR basophil OR allergy OR asthma)',
 "NK/T cell":
   'lipoteichoic acid AND ("NK cell" OR "T cell" OR splenocyte OR Th1 OR Th17)',
 "epithelial/lung":
   'lipoteichoic acid AND (epithelial OR alveolar OR lung OR airway)',
 "neutrophil":
   'lipoteichoic acid AND (neutrophil OR granulocyte OR "whole blood")',
 # structural chemistry
 "D-alanine content":
   'lipoteichoic acid AND ("D-alanine" OR "D-alanyl" OR alanylation)',
 "glycolipid anchor":
   'lipoteichoic acid AND (glycolipid OR "anchor" OR diglucosyl OR Glc2DAG)',
 "chain length":
   'lipoteichoic acid AND ("chain length" OR "repeating unit" OR polymerization)',
 "fatty acid composition":
   'lipoteichoic acid AND ("fatty acid" OR acyl OR "lipid moiety" OR deacylation)',
 # genomics-era vocabulary (2006-)
 "WGS/comparative genomics":
   '(lipoteichoic OR "teichoic acid") AND ("whole genome" OR "comparative genomics" '
   'OR "genome sequence" OR "genome mining")',
 "in silico prediction":
   '("teichoic acid" OR lipoteichoic) AND ("in silico" OR bioinformatic* OR "gene cluster" '
   'OR operon) AND (predict* OR annotat*)',
 "pangenome":
   '(Lactobacillus OR Lactiplantibacillus OR Streptococcus) AND (pangenome OR "pan-genome") '
   'AND ("cell wall" OR teichoic)',
 "dlt operon genomics":
   '"dlt operon" AND (genome OR sequence OR distribution OR phylogen*)',
 "probiotic genomics":
   'probiotic AND ("genome-based" OR "genomic analysis" OR "strain selection") '
   'AND (immunomodulat* OR "cell wall")',
 # other readouts
 "antitumour":
   'lipoteichoic acid AND (antitumor OR antitumour OR tumor necrosis)',
 "biofilm/adhesion":
   'lipoteichoic acid AND (biofilm OR adhesion OR autolysis)',
 "sepsis/shock":
   'lipoteichoic acid AND (sepsis OR "septic shock" OR endotoxemia)',
 "commercial/contamination":
   'lipoteichoic acid AND (commercial OR contaminat* OR "polymyxin B" OR endotoxin-free '
   'OR lipopeptide)',
}

# ---- inclusion criteria applied to every abstract -----------------------
CRITERIA = [
    "the paper reports purified lipoteichoic acid, not whole cells",
    "at least two strains are compared",
    "the comparison is within one assay system",
    "an activity grade can be read for each strain",
]


ALL_QUERIES = {**QUERIES_METHOD, **QUERIES}   # 14 + 25 = 39


def run():
    out = {}
    for name, q in ALL_QUERIES.items():
        for attempt in range(5):
            try:
                d = json.loads(urllib.request.urlopen(
                    ESEARCH + urllib.parse.quote(q), timeout=60).read().decode())
                break
            except Exception:
                time.sleep(2 + 2 * attempt)
        out[name] = {"query": q,
                     "total": int(d["esearchresult"]["count"]),
                     "ids": d["esearchresult"]["idlist"]}
        print(f"{name:26s} total={out[name]['total']:>6d}")
        time.sleep(0.45)
    json.dump(out, open("query_results.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    run()
