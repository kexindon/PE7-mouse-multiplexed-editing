# PE7-mouse-multiplexed-editing
Analysis associated with PE7 manuscript.

## Figure legends (Updated Sep 24, 2026):

**PE7 multiplexed prime editing recapitulates the top-VAF driver mutation landscape in the majority of patients.**

Patient-level coverage of the highest-VAF somatic mutations across the AACR Project GENIE cohort (v15.0). For each tumor sample, mutations were ranked by variant allele frequency and the top three retained. A patient is scored as "covered" at rank N only if all of their top-1 to top-N mutations are individually amenable to the indicated modeling strategy; denominators therefore comprise the patients carrying at least N VAF-annotated mutations (Single, n = 131,529; Dual, n = 112,425; Triple, n = 96,170). Gray, all eligible patients; blue, mutations with a precisely modelable mouse orthologous allele; purple, mutations addressable by any breeding existing genetically engineered mouse models (GEMMs); orange, alleles for which at least one NGG-PAM pegRNA could be designed by PEGG for the mouse allele. Point = an exact gene-plus-amino-acid match to a published point-mutation/knock-in/humanized allele; LoF = a knockout or conditional allele exists for a truncating loss-of-function mutation.LoF mutations include Nonsense, Frameshift, Splice_Site, Splice_Region, Nonstop and Translation_Start_Site variant classes. VAF = Variant Allele Frequency.

## Methods (Updated Sep 24, 2026): 

**Estimation of patient-level modelability of top-VAF driver combinations**

Patient mutations and clinical annotations were obtained from AACR Project GENIE v15.0 (https://genie.synapse.org/Explore/GENIE). Gene symbols were harmonized to current HGNC nomenclature prior to joining. Each human mutation was matched to mouse mutation by using H2M (REF: https://www.nature.com/articles/s41587-025-02925-0). For every H2M-modelable mouse allele, pegRNAs were designed with PEGG (REF: https://www.nature.com/articles/s41587-024-02172-9) against the GRCm39 mouse reference assembly (GCF_000001635.27), using an NGG PAM, a primer binding site length of 7 nt and a reverse transcriptase template length of 15 nt.

Coverage by existing genetically engineered mouse models. We looked at two categories. First, the mutation matched an existing point-mutation, knock-in, humanized, or spontaneous/chemically induced allele at both the gene and the precise amino-acid substitution. Second, a knockout or conditional allele existed for the geneand the patient mutation was a truncating loss-of-function variant (Nonsense_Mutation, Frame_Shift_Del, Frame_Shift_Ins, Splice_Site, Splice_Region, Nonstop_Mutation or Translation_Start_Site). We assume that any combination of existing alleles can in principle be bred into a single animal.
