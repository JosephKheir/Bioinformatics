Overview
--------

16S metagenomic next generation sequencing (mNGS) is a targeted method
in microbial genomics for pathogen and microbial/ microbiome analysis
and profiling. The general approach uses primers that anneal to
universally conserved regions on the 16S sequence (Chiu and Miller
2019). Through PCR, these DNA sequences are amplified and then used for
downstream analysis, such as in Qiime2 mNGS micribiome analysis (Bolyen
et al. 2019). In this module, QIIME2 workflow piplines were analyzed
using two mNGS studies and their datasets as demonstration.

Methods
-------

### Fecal Microbiota Transplat (FMT) Study

Autism spectrum disorder (ASD) is a neurobiological disorder evident by
characteristicly impaired communication and social interactions. ASD is
generally characterized by restricted, repetitive stereotyped patters of
behavior and interests (Kang et al. 2017). In this tutorail, the data
used was a subset of the small open label study where the impact of
Microbiota Transfer Therapy (MTT) on gut microbiota was evaluated in 18
ASD diagnosed children between the ages of 7 - 17 years old. Gut
microbiota composition, GI symptoms and ASD clinical symptoms were
evaluated. Each participant participated in the study for a duration of
18 weeks with a 10 week treatment period and an 8 week follow-up period.
MTT was administered either orally or rectally. 20 age and gender
matched controls were selected for the study as well. Microbial DNA was
extracted and isolated from subjects and DNA amplification was performed
using next generation sequencing mNGS. Sample metadata was collected
capturing datapoints such as sample ids, sample type, treatment group,
subject id, sample collection week, route of administration, subject
age, gender, height and weight. Barcode data was also collected for each
sequencing run. Sample metadata was processed using Quantitative
Insights Into Microbial Ecology (QIIME) and microbial population and
taxonomic characteristics were analyzed (Kang et al. 2017).

### Atacama Soil Microbiome Study

Desert regions are globally responsible for approximately 27% of soil
organic carbon storage (<span class="citeproc-not-found"
data-reference-id="Nielson">**???**</span>). The continued degradation
due to climate change is contributing to a loss in productivity in these
vital ecosystems. This degradation contributes to disruptions in
microbial communities of desert microbiomes lending to decreased
phylogenetic diversity and functions (eg. nutrient cycling). In this
data set, soil samples from were collected in two parallel (east - west)
transects of the Atacama Desert of norther Chile. The transects spanned
through both hyperarid and arid regions of the desert. The northern most
parallel transect was referred to as Baquedano and the southern most
parallel transect was referred to as Yungay. 10 and 12 sites were
located along the norther and souther transects respectively. Each
collection site was assigned a site name/ id. In addition to soil
samples, at each site, metadata points such as elevation, relative soil
humidity (RH), collection site soil ph, soil temperature and vegetation
cover were collected. Each collection was labeled with a sample id and
assigned a sequence barcode for identification. Additionally, sample
extract concentration as well as amplicon concentration was recorded per
collection sample. Qiime bioinformatics analysis tools were employed to
generate beta diversity profiles of the Atacama desert relative to the
cross- biome survey of global soil samples (Neilson et al. 2017).

References
----------

Bolyen, Evan, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich,
Christian C. Abnet, Gabriel A. Al-Ghalith, Harriet Alexander, et al.
2019. “Reproducible, Interactive, Scalable and Extensible Microbiome
Data Science Using QIIME 2.” *Nature Biotechnology* 37 (8): 852–57.
<https://doi.org/10.1038/s41587-019-0209-9>.

Chiu, Charles Y., and Steven A. Miller. 2019. “Clinical Metagenomics.”
*Nature Reviews Genetics* 20 (6): 341–55.
<https://doi.org/10.1038/s41576-019-0113-7>.

Kang, Dae-Wook, James B. Adams, Ann C. Gregory, Thomas Borody, Lauren
Chittick, Alessio Fasano, Alexander Khoruts, et al. 2017. “Microbiota
Transfer Therapy Alters Gut Ecosystem and Improves Gastrointestinal and
Autism Symptoms: An Open-Label Study.” *Microbiome* 5 (1): 10.
<https://doi.org/10.1186/s40168-016-0225-7>.

Neilson, Julia W., Katy Califf, Cesar Cardona, Audrey Copeland, Will van
Treuren, Karen L. Josephson, Rob Knight, et al. 2017. “Significant
Impacts of Increasing Aridity on the Arid Soil Microbiome.” *mSystems* 2
(3): e00195–16. <https://doi.org/10.1128/mSystems.00195-16>.
