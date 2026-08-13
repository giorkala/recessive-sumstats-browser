# Frequently Asked Questions

### What is BRaVa?

The **Biobank Rare Variant Analysis (BRaVa) consortium** harmonises rare coding-variant association analyses across global biobanks and meta-analyses them, maximising power and ancestral diversity for gene-based rare-variant discovery.

### What study are these summary statistics from?

These results are from the BRaVa study **“Meta-analysis across six global biobanks identifies recessive coding associations with complex traits and diseases.”**

The study analysed recessive gene-based effects across up to **948,690 exome- or whole-genome-sequenced participants** and 41 phenotypes.

### Which cohorts are included?

The analysis includes six biobanks:

- **AOU** — All of Us Research Program
- **BioMe** — BioMe Biobank
- **BBJ** — BioBank Japan
- **G&H / GNH** — Genes & Health
- **UKB** — UK Biobank
- **100kGP / GEL** — Genomics England 100,000 Genomes Project

Analyses were generally carried out separately within genetically inferred ancestry groups before meta-analysis.

### What does “recessive” mean in these analyses?

For each gene, individuals were classified according to whether qualifying variants affected **zero, one, or both copies of the gene**.

The recessive model compares individuals with a qualifying **bi-allelic genotype** against everyone else. In other words, the recessive genotype is encoded as:

`0, 0, 1`

for zero, one, or two affected gene copies respectively.

### What is a bi-allelic genotype?

A bi-allelic genotype means that **both copies of a gene contain qualifying variants**.

This can occur through:

- **Homozygosity:** the same qualifying variant is present on both copies of the gene.
- **Compound heterozygosity (CH):** two different qualifying variants occur on opposite copies of the gene.

Statistical phasing was used to determine whether different variants were on opposite haplotypes and could therefore form a compound-heterozygous genotype.

### Were compound heterozygous genotypes included?

Yes, where phased data were available. Including compound heterozygotes increased the number of detectable damaging bi-allelic genotypes by approximately **19%**.

**All of Us was not phased for this analysis**, so compound-heterozygous genotypes were not identified in AOU; its recessive analyses therefore rely on homozygous genotypes.

### What MAF threshold was used?

Qualifying variants were restricted to variants with **minor allele frequency (MAF) <5%**.

### What do the annotation masks mean?

Four gene-level variant masks were considered.

**pLoF**

Predicted loss-of-function variants. These represent variants predicted to severely disrupt gene function.

**pLoF | damaging_missense**

Includes pLoF variants plus predicted damaging missense/protein-altering variants.

Damaging variants include low-confidence loss-of-function variants and variants meeting one or more of the following prediction thresholds:

- **REVEL ≥ 0.773**
- **CADD ≥ 28.1**
- splice variants with **SpliceAI Δ ≥ 0.5**



**nonsynonymous**

The broadest protein-altering mask. It includes:

- pLoF variants
- damaging missense/protein-altering variants
- other missense/protein-altering variants



**synonymous**

Synonymous coding variants. These generally do not alter the encoded amino acid and were included primarily as a **negative-control mask**.

### Are the annotation masks nested?

Yes. Conceptually:

`pLoF ⊂ pLoF|damaging_missense ⊂ nonsynonymous`

The broader masks therefore contain progressively larger sets of qualifying variants.

### What is the additive model shown in the browser?

The corresponding additive analysis encodes the number of **affected gene copies** as:

`0, 1, 2`

This is slightly different from a conventional rare-variant additive burden based on the total number or weighted sum of qualifying variants carried by an individual.

### How were ancestry groups handled?

Each biobank was divided into genetically inferred ancestry groups. Gene–phenotype association testing was performed within these **biobank × ancestry subcohorts**, which were subsequently combined in meta-analysis.

The ancestry groups represented across the study include **AFR, AMR, EAS, EUR and SAS**, although not every ancestry is represented in every biobank.

### How many individuals with bi-allelic genotypes were required for a test?

An association was evaluated only when at least **five individuals with a bi-allelic genotype** were available for that gene, phenotype and subcohort.

### How were results meta-analysed?

Association testing was first performed within individual biobank–ancestry subcohorts.

For each damaging annotation mask, evidence was then combined across subcohorts using meta-analysis. The study also combined evidence across the different nonsynonymous masks using the **Cauchy Combination Test (CCT)** to obtain an overall gene–trait association p-value.

### What does the effect estimate mean?

For **binary phenotypes**, effect sizes are reported as odds ratios or their corresponding regression estimates.

For **quantitative phenotypes**, effect estimates represent changes in the transformed phenotype, generally expressed in standard-deviation units.

The direction of effect refers to individuals carrying a bi-allelic qualifying genotype compared with individuals without one.

### Why can effect estimates be very large?

Recessive qualifying genotypes can be extremely rare. Consequently, some gene–phenotype tests are based on relatively few bi-allelic carriers.

Large effect estimates and wide confidence intervals should therefore be interpreted together with:

- the number of carriers,
- the contributing cohorts,
- consistency of effects between cohorts, and
- the statistical uncertainty.

### Why are some cohort-level results missing?

A gene–phenotype combination may not have been tested in every cohort or ancestry. Possible reasons include:

- insufficient numbers of bi-allelic carriers,
- the phenotype not being available in that cohort,
- ancestry-specific sample availability, or
- cohort-specific QC or analysis requirements.

A missing result should therefore **not be interpreted as evidence of no association**.

### Why do some genes have results for one annotation but not another?

Each annotation mask contains a different set of variants. A gene may therefore have enough bi-allelic carriers to test a broad mask such as `nonsynonymous`, but too few carriers to test a stricter mask such as `pLoF`.

### What is the difference between a cohort-level result and a meta-analysis result?

**Cohort-level results** show the association estimated within an individual biobank and ancestry group.

**Meta-analysis results** combine evidence across contributing subcohorts, giving greater statistical power and allowing evidence to be assessed across populations.

### Does a significant association prove a recessive mechanism?

No. A significant recessive association indicates that individuals with bi-allelic qualifying genotypes differ in the phenotype.

The study additionally compared recessive and additive models to identify associations whose evidence was substantially stronger under a recessive model. The paper classified an association as more consistent with recessive inheritance when the recessive p-value was more than 100-fold smaller than the corresponding additive p-value.

### How should these results be interpreted?

The browser is intended as a resource for exploring gene-based recessive associations.

Associations should be interpreted alongside the number of carriers, effect estimates, confidence intervals, consistency across cohorts and ancestries, annotation mask, biological evidence, and results reported in the publication.

These summary statistics should not by themselves be interpreted as evidence that a particular variant is pathogenic or clinically actionable.

### Where can I find the full methodology?

Full details of cohort processing, phenotype definitions, variant annotation, phasing, association testing and meta-analysis are available in the accompanying BRaVa publication and supplementary information.

The browser provides a convenient interface to the released cohort-level and meta-analysis summary statistics.