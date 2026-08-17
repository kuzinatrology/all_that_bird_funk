<h1 align="center">all_that_bird_funk</h1>

<p align="center">
  <b>Comparative genomics of avian longevity</b><br>
  Signatures of selection in ageing-related genes across 107 bird species
</p>

---

Birds live two to three times longer than mammals of comparable body mass. Long life has arisen repeatedly in unrelated
avian lineages. This pipeline searches for episodic diversifying selection in
ageing-associated genes across those lineages and relates it to size-corrected
longevity.

Maximum lifespan scales with body mass, so the analysis uses LQ. An allometric
regression log₁₀(L) = β₀ + β₁·log₁₀(M) is fitted across the sampled species, and LQ is
observed lifespan divided by predicted lifespan. The top 20% of LQ form the
**foreground** branch set for the branch-site tests.

## Pipeline

```mermaid
flowchart TD
    A["Bird proteomes<br/>107 species, longest isoform per gene"] --> B["Homology search<br/>MMseqs2 easy-search"]
    Q["Query proteins<br/>OpenGenes · Matrisome · regeneration"] --> B
    B --> C["CDS retrieval<br/>one FASTA per gene"]
    C --> D["Power gate<br/>≥18 species, ≥4 long-lived"]
    D --> E["Codon alignment<br/>MACSE alignSequences"]
    E --> F["Frameshifts and stops → gap codons<br/>MACSE exportAlignment"]
    F --> G["Mask segments off the alignment profile<br/>HmmCleaner"]
    G --> H["Project mask onto codons<br/>MACSE reportMaskAA2NT"]
    H --> I["Coverage QC<br/>realign affected genes"]
    I -.->|realign| E
    H --> J["Drop sequences that collapsed under masking,<br/>then sparse codon columns"]
    J --> K["Prune and label species tree<br/>ete3"]
    K --> L["BUSTED-E on foreground branches<br/>all genes"]
    L --> N["Benjamini-Hochberg<br/>q &lt; 0.2"]
    N --> O["BUSTED-E on background branches"]
    O --> P["No signal on the background:<br/>selection confined to long-lived"]
    P --> R["RELAX: K on the foreground"]
    P --> S["aBSREL: which branches"]
    P --> T["MEME: which codons"]
    P --> U["FitMG94: per-species ω"]
    U --> M["PGLS: LQ on ω"]
```

## Input data

| | |
|---|---|
| **Sequences** | proteomes, CDS and genome annotations for the bird species |
| **Traits** | body mass and maximum lifespan per species, from which LQ is derived |
| **Queries** | amino acid sequences of genes of interest: OpenGenes ageing categories, Matrisome Project, tissue regeneration sets |
| **Phylogeny** | a species tree covering the sampled species |

## Steps

| # | Step | Tool | Notes |
|---|---|---|---|
| 1 | Homology search | MMseqs2 `easy-search` | `-c 0.2 --min-seq-id 0.3 --alt-ali 25`, best hit per species. The target database keeps one longest isoform per gene, so a gene is represented once per species |
| 2 | CDS retrieval | — | the coding sequence behind each protein hit, collected per gene |
| 3 | Cleaning and gate | — | species names reduced to binomials; fragments shorter than 10% of the gene mean dropped; a gene is kept only with ≥18 species and ≥4 long-lived species |
| 4 | Codon alignment | MACSE `alignSequences` | frameshift-aware alignment of coding sequences |
| 5 | Export | MACSE `exportAlignment` | frameshifts and stop codons, including the terminal one, become whole gap codons, since HyPhy rejects alignments containing stop codons |
| 6 | Masking | HmmCleaner `--large --specificity` | scores every residue against an HMM profile built from the alignment itself and masks segments that do not fit it. In practice these are mispredicted exons and frameshifted stretches |
| 7 | Mask transfer | MACSE `reportMaskAA2NT` | the amino-acid mask is projected onto the nucleotide alignment |
| 8 | Coverage QC | — | a sequence delivering under 30% of the codons of its own CDS is removed from the inputs and its gene is realigned |
| 9 | Filtering | — | a sequence retaining under half of its own length after masking is dropped from that gene. Then a codon column occupied in under 50% of the remaining sequences is deleted |
| 10 | Trees | ete3 | the species tree is pruned to each gene's species and labelled `{Foreground}` / `{Background}` |
| 11 | Selection, all genes | HyPhy `busted` | BUSTED-E on the foreground branches, with an error sink absorbing alignment artefacts. P-values corrected with Benjamini-Hochberg; genes with q < 0.2 go on |
| 12 | Confinement to long-lived lineages | HyPhy `busted` | the same genes are rerun with the long-lived species as background. A gene with no background signal carries selection confined to the long-lived lineages, and only those genes go on |
| 13 | Selection, confined genes | HyPhy | RELAX: is selection intensified on the foreground (K). aBSREL: which branches carry it. MEME: which codons, and which lineages drive each. FitMG94 `--type lineage`: per-species root-to-tip ω |
| 14 | Association | PGLS | LQ regressed on per-species ω, Brownian-motion covariance from the pruned tree |

The cleaning scheme follows Botero-Castro et al. (2017), with two changed thresholds:
column occupancy 50% (60% in the paper) and HmmCleaner in `--specificity` mode. The
HmmCleaner authors report specificity dropping in fast-evolving regions, which is
where the signal here is. Measured on this dataset, that mode masks 0.8% of residues
and shifts the column filter outcome by 0.04%.

Every step that deletes columns records an old→new index map, so site coordinates
reported by MEME can be translated back to the original alignment.

## Running the selection analyses

Every HyPhy call passes `ENV=TOLERATE_NUMERICAL_ERRORS=1`. Without it a fraction of
genes aborts with `Internal error in ComputeBranchCache`, which HyPhy reports as
numerical instability on large alignments.

```bash
hyphy CPU=1 ENV=TOLERATE_NUMERICAL_ERRORS=1; busted \
    --alignment <gene>_binomial_NT.fasta \
    --tree      <gene>_labeled.nwk \
    --branches  Foreground \
    --error-sink Yes \
    --output    <gene>_busted.json
```

`busted_e_run.py` runs this across the whole gene set: one process per core, a worker
cap adjustable while the run is in flight, a per-gene timeout, and back-off when other
users load the machine.

## Limitations

- Bird sequences are best MMseqs2 hits per human query. Orthology was not inferred
  formally.
- 342 of 3403 genes are excluded as too long to align (median CDS 5313 nt against
  1275 for the rest; MACF1 is 22 479 nt). The gene set is biased against long CDS.
- DIO2, SELENOP and MT-CYB are run separately. Selenocysteine is encoded by TGA,
  which MACSE reads as a stop codon; MT-CYB needs `-gc_def 2`.

## Contents

| File | Role |
|---|---|
| `pipeline.ipynb` | the main notebook: all stages, from homology search to final tables |
| `cicl_mmseq_easy_search.sh` | homology search |
| `build_primary_merged.py` | builds the one-isoform-per-gene target database |
| `qc_redo_lowcov.py` | coverage QC and realignment of affected genes |
| `rebuild_labeled_trees.py` | pruning and labelling of per-gene trees |
| `busted_e_run.py` | BUSTED-E across the full gene set |
| `relax_run.sh`, `absrel.sh`, `run_meme.sh`, `run_fitmg94.sh`, `run_busted_background.sh` | the remaining HyPhy analyses |
| `pruned_tree.nwk` | species tree |

## References

- Ranwez V., Douzery E.J.P., Cambon C., Chantret N., Delsuc F. **MACSE v2.** *Mol Biol Evol* 2020.
- Di Franco A., Poujol R., Baurain D., Philippe H. **HmmCleaner.** *BMC Evol Biol* 2019.
- Botero-Castro F., Figuet E., Tilak M.-K., Nabholz B., Galtier N. **Avian genomes revisited.** *Mol Biol Evol* 2017.
- Kosakovsky Pond S.L. et al. **HyPhy 2.5.** *Mol Biol Evol* 2020. · Murrell B. et al. **BUSTED.** *Mol Biol Evol* 2015.
- Steinegger M., Söding J. **MMseqs2.** *Nat Biotechnol* 2017.
