# HypoHyperMDD
Differential gene expression analysis of hyperphagic and hypophagic MDD

**Abstract**

Major Depressive Disorder (MDD) is a heterogenous and etiologically complex disease encompassing a broad spectrum of psychopathology, presumably arising from distinct pathophysiological mechanisms. Divergent appetitive phenotypes including Hyperphagic MDD (characterized by an increased appetite) and Hypophagic MDD (characterized by a decrease in appetite) are important clinical characteristics that are closely related to comorbidities, including cardiometabolic disorders. Prior evidence supports the notion that hyperphagia is associated with atypical depression, decreased stress-hormone signaling, a pro-inflammatory status, hypersomnia, and poorer clinical outcomes. Yet, our understanding of the underlying mechanisms of Hyperphagic and Hypophagic MDD is limited, and knowledge of associated biological correlates of these endophenotypes remain fragmented. We performed an exploratory study on peripheral blood RNA profiling using bulk RNAseq in unmedicated individuals with Hyperphagic and Hypophagic MDD (n=8 and n=13, respectively) and discovered individual genes and gene pathways associated with appetitive phenotypes. In addition, we used the Maastricht Acute Stress Task to uncover stress-related transcriptomic profiles in Hyper- and Hypophagic MDD. 

**Data and Analysis**
- **Raw Data**: RNAseq fastq files can be accessed through GEO accession GSE231347 (token protected as of Sept 17th, 2023)
- **Metadata**:
  - de_Kluiver_ref_pathways_all.csv: Reference pathway mentioned in de Kluiver et al., used for GSEA
  - metadata_v2_goi_count_t0.csv and metadata_v2_goi_count_t105.csv: deidentified metadata for RNAseq
- **Analysis**: Main analysis was done in Analysis.Rmd
