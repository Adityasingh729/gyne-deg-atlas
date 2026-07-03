# Dataset Summary & Mapping Report

This report compiles the filled screening and mapping metadata for the analyzed datasets in the Gyne-DEG Atlas project across PCOS, Endometriosis, and Adenomyosis.

## Master Dataset Table

| A1 GEO ID | A2 Disease+tissue | A3 GEO title | B1 Case group | B2 Control group | B3 Control type | B4 n vs n | C1 Type | C2 Arm used | D1 Data | D2 ID type | D3 SRP | E1 Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| GSE207000 | Adenomyosis / endometrium | Transcriptome analysis of EGR1 knock-down eutopic endometrial mesenchymal stem cells by RNA-sequencing in patients with adenomyosis. | EGR1-KD | Scr | scrambled shRNA control | 3 v 3 | RNA-seq | EGR1-KD vs Scr | raw counts | Gene Symbol | SRP383569 | ANALYSED - DESeq2; EGR1 knockdown in eutopic MSCs, strong knockdown signature (5049 DEGs) |
| GSE262302 | Adenomyosis / plasma exosomes | Exosomal miR-92a-3p as a diagnostic biomarker for adenomyosis and its pathological mechanism | Adenomyosis | Control | healthy control | 8 v 8 | RNA-seq | Adenomyosis vs Control | raw counts | Gene Symbol | SRP498967 | ANALYSED - DESeq2; plasma exosomes miRNA profile, low signal (7 DEGs) |
| GSE316028 | Endometriosis / endometrial stromal cells | Pyruvate carboxylase promotes glycolysis and progression of endometriosis by activating the AKT pathway | siPC | NC | negative control siRNA | 3 v 3 | RNA-seq | siPC vs NC | raw counts | Gene Symbol | SRP568341 | ANALYSED - DESeq2; PC knockdown in immortalized ESCs, very strong signature (4643 DEGs) |
| GSE40007 | Endometriosis / endometrial stromal cells | TNFalpha and IL1beta stimulate differential gene expression in endometrial stromal cells | TNFα | Vehicle | vehicle control | 3 v 3 | Microarray | TNFα vs Vehicle | microarray intensities | Entrez | no srp | ANALYSED - limma; TNFa stimulated ESCs, low response (2 DEGs) |
| GSE44207 | Endometriosis / endometrial stromal cells | CCAAT/enhancer-binding protein alpha is epigenetically silenced by histone acetylation in endometriosis and promotes the pathogenesis of endometriosis: a novel therapeutic target | VPA-treated | Untreated | untreated control | 4 v 4 | Microarray | VPA-treated vs Untreated | microarray intensities | Probe ID | no srp | ANALYSED - limma; VPA treatment of primary ESCs, massive chromatin remodeling response (3895 DEGs) |
| GSE47360 | Endometriosis / endometrial stromal cells | Genome-wide DNA methylation profiling in cultured eutopic and ectopic endometrial stromal cells (expression) | Endometriosis | Control | healthy control | 6 v 3 | Microarray | Endometriosis vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; ESC case vs control, no independent signal (flag weak) |
| GSE56854 | Endometriosis / endometrial stromal cells | Compulsory expression of miR-210 in normal endometrial stromal cells | miR-210 | NC miRNA | negative control transfection | 4 v 4 | Microarray | miR-210 vs NC miRNA | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; miR-210 overexpression in ESCs, strong response (1130 DEGs) |
| GSE75425 | Endometriosis / endometrial stromal cells | mRNA expression profiles in decidualized and non-decidualized normal endometrial stromal cells (NESCs). | cAMP + Dienogest | FBS | vehicle control | 4 v 4 | Microarray | cAMP + Dienogest vs FBS | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; decidualized ESC cAMP/dienogest vs FBS, moderate signature (6 DEGs) |
| GSE51981 | Endometriosis / endometrium | Molecular Classification of Endometriosis and Disease Stage Using High-Dimensional Genomic Data | Endometriosis | Control | healthy control | 77 v 71 | Microarray | Endometriosis vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; bulk endometrium case vs control, massive diagnostic transcriptomic signature (8476 DEGs) |
| GSE6364 | Endometriosis / endometrium (eutopic) | Endometriosis eutopic endometrium microarray expression | Endometriosis | Control | healthy control | 21 v 16 | Microarray | Endometriosis vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; secretory phase endometrium, zero sig genes (flag weak) |
| GSE286315 | PCOS / endometrial stromal cells | HOXC4 Promotes Proliferation of Endometriotic Stromal Cells via the SLIT2-ROBO1 Axis | HOXC4 KO | Control EcSCs | non-targeting control | 3 v 3 | RNA-seq | HOXC4 KO vs Control EcSCs | raw counts | Gene Symbol | SRP518602 | ANALYSED - DESeq2; HOXC4 knockdown in endometriotic stromal cells, moderate signature (82 DEGs) |
| GSE102293 | PCOS / granulosa | Enhanced Inflammatory Transcriptome in the Granulosa Cells of Women With Polycystic Ovarian Syndrome | PCOS | Control | healthy control | 2 v 4 | Microarray | PCOS vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; granulosa cells, small sample size (2 v 4) yields zero sig (flag weak) |
| GSE106724 | PCOS / granulosa | Profiles for long noncoding RNAs in ovarian granulosa cells from polycystic ovary syndrome patients with different serum concentrations of androgen | PCOS (HA) | Control | healthy control | 4 v 4 | Microarray | PCOS (HA) vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; PCOS hyperandrogenism in granulosa, strong signature (1088 DEGs) |
| GSE114419 | PCOS / granulosa | Coding and Non-coding gene expression signatures of granulosa cells from Polycystic Ovary Syndrome Patients | PCOS | Control | healthy control | 3 v 3 | Microarray | PCOS vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; granulosa cells, small sample size (3 v 3), zero sig (flag weak) |
| GSE137684 | PCOS / granulosa | Differential messenger RNA expression in Granulosa Cells from polycystic ovary syndrome with Normoandrogen and Hyperandrogen: Identification of gene sets through bioinformatic Filtering analysis | PCOS (NA) | Control | healthy control | 4 v 4 | Microarray | PCOS (NA) vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; normoandrogenic PCOS in granulosa, zero sig (flag weak) |
| GSE34526 | PCOS / granulosa (whole-gc) | Differential Gene Expression in Granulosa Cells from Polycystic Ovary Syndrome Patients with and without Insulin Resistance: Identification of Susceptibility Gene Sets through Network Analysis | PCOS | Control | healthy control | 7 v 3 | Microarray | PCOS vs Control | microarray intensities | Probe ID | no srp | ANALYSED - limma; granulosa cells PCOS vs control, zero sig (flag weak) |
| GSE6798 | PCOS / skeletal muscle | Reduced expression of mitochondrial oxidative metabolism genes in skeletal muscle of women with PCOS | PCOS | Control | healthy control | 16 v 13 | Microarray | PCOS vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; skeletal muscle, low signal (21 DEGs) |
| GSE8157 | PCOS / skeletal muscle | Gene expression profiling in skeletal muscle of PCOS after pioglitazone therapy | PCOS | Control | healthy control | 10 v 13 | Microarray | PCOS vs Control | microarray intensities | Gene Symbol | no srp | ANALYSED - limma; skeletal muscle PCOS vs control, strong signature (1723 DEGs) |
| GSE199225 | PCOS / skeletal muscle (vastus) | TGFβ1 impairs the transcriptomic response to contraction in myotubes from women with polycystic ovary syndrome | PCOS | Control | healthy control | 24 v 25 | RNA-seq | PCOS vs Control | raw counts (recount3) | Gene Symbol | SRP365287 | ANALYSED - DESeq2; vastus skeletal muscle, robust differential response (995 DEGs) |

---

## Individual Dataset Sample Sheets (Sample ID -> Group -> Tissue)

### Sample Sheet for GSE207000

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM6267779 | Case | EGR1-KD | Eutopic endometrium from patients with adenomyosis |
| GSM6267780 | Case | EGR1-KD | Eutopic endometrium from patients with adenomyosis |
| GSM6267781 | Case | EGR1-KD | Eutopic endometrium from patients with adenomyosis |
| GSM6267782 | Control | Scr | Eutopic endometrium from patients with adenomyosis |
| GSM6267783 | Control | Scr | Eutopic endometrium from patients with adenomyosis |
| GSM6267784 | Control | Scr | Eutopic endometrium from patients with adenomyosis |

### Sample Sheet for GSE262302

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM8163035 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163036 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163037 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163038 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163039 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163040 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163041 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163042 | Case | Adenomyosis | adenomyosis patients plasma exosomes |
| GSM8163043 | Control | Control | non-adenomyosis patients plasma exosomes |
| GSM8163044 | Control | Control | non-adenomyosis patients plasma exosomes |
| GSM8163045 | Control | Control | non-adenomyosis patients plasma exosomes |
| GSM8163046 | Control | Control | non-adenomyosis patients plasma exosomes |
| GSM8163047 | Control | Control | non-adenomyosis patients plasma exosomes |
| GSM8163048 | Control | Control | non-adenomyosis patients plasma exosomes |
| GSM8163049 | Control | Control | non-adenomyosis patients plasma exosomes |
| GSM8163050 | Control | Control | non-adenomyosis patients plasma exosomes |

### Sample Sheet for GSE316028

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM9443541 | Case | siPC | immortalized human endometriotic stromal cells |
| GSM9443542 | Case | siPC | immortalized human endometriotic stromal cells |
| GSM9443543 | Case | siPC | immortalized human endometriotic stromal cells |
| GSM9443538 | Control | NC | immortalized human endometriotic stromal cells |
| GSM9443539 | Control | NC | immortalized human endometriotic stromal cells |
| GSM9443540 | Control | NC | immortalized human endometriotic stromal cells |

### Sample Sheet for GSE40007

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM982999 | Control | vehicle | human endometrial stromal cell line (T-HESC) (ATCC CRL-4003) |
| GSM983000 | Control | vehicle | human endometrial stromal cell line (T-HESC) (ATCC CRL-4003) |
| GSM983001 | Control | vehicle | human endometrial stromal cell line (T-HESC) (ATCC CRL-4003) |
| GSM983002 | Case | TNFα | human endometrial stromal cell line (T-HESC) (ATCC CRL-4003) |
| GSM983003 | Case | TNFα | human endometrial stromal cell line (T-HESC) (ATCC CRL-4003) |
| GSM983004 | Case | TNFα | human endometrial stromal cell line (T-HESC) (ATCC CRL-4003) |

### Sample Sheet for GSE44207

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM1080771 | Control | untreated | ECSC, untreated, replicate 1 |
| GSM1080772 | Control | untreated | ECSC, untreated, replicate 2 |
| GSM1080773 | Control | untreated | ECSC, untreated, replicate 3 |
| GSM1080774 | Control | untreated | ECSC, untreated, replicate 4 |
| GSM1080775 | Case | VPA-treated | ECSC, VPA-treated, replicate 1 |
| GSM1080776 | Case | VPA-treated | ECSC, VPA-treated, replicate 2 |
| GSM1080777 | Case | VPA-treated | ECSC, VPA-treated, replicate 3 |
| GSM1080778 | Case | VPA-treated | ECSC, VPA-treated, replicate 4 |

### Sample Sheet for GSE47360

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM1148101 | Control | non-endometriosis | eutopic endometrium |
| GSM1148102 | Control | non-endometriosis | eutopic endometrium |
| GSM1148103 | Control | non-endometriosis | eutopic endometrium |
| GSM1148104 | Case | endometriosis | eutopic endometrium |
| GSM1148105 | Case | endometriosis | eutopic endometrium |
| GSM1148106 | Case | endometriosis | eutopic endometrium |
| GSM1148107 | Case | endometriosis | ovarian chocolate cyst |
| GSM1148108 | Case | endometriosis | ovarian chocolate cyst |
| GSM1148109 | Case | endometriosis | ovarian chocolate cyst |

### Sample Sheet for GSE56854

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM1370215 | Control | normal endometrial stromal cells transfected with negative control precursor miRNA | Human proliferative-phase endometrium |
| GSM1370216 | Case | normal endometrial stromal cells transfected with precursor hsa-miR-210 | Human proliferative-phase endometrium |
| GSM1370217 | Control | normal endometrial stromal cells transfected with negative control precursor miRNA | Human proliferative-phase endometrium |
| GSM1370218 | Case | normal endometrial stromal cells transfected with precursor hsa-miR-210 | Human proliferative-phase endometrium |
| GSM1370219 | Control | normal endometrial stromal cells transfected with negative control precursor miRNA | Human proliferative-phase endometrium |
| GSM1370220 | Case | normal endometrial stromal cells transfected with precursor hsa-miR-210 | Human proliferative-phase endometrium |
| GSM1370221 | Control | normal endometrial stromal cells transfected with negative control precursor miRNA | Human proliferative-phase endometrium |
| GSM1370222 | Case | normal endometrial stromal cells transfected with precursor hsa-miR-210 | Human proliferative-phase endometrium |

### Sample Sheet for GSE75425

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM1954914 | Control | treatment: 12d 10% charcoal-stripped heat-inactivated FBS | Proliferative phase normal endometrium, 12 days, replicate1 |
| GSM1954915 | Control | treatment: 12d 10% charcoal-stripped heat-inactivated FBS | Proliferative phase normal endometrium, 12 days, replicate2 |
| GSM1954916 | Control | treatment: 12d 10% charcoal-stripped heat-inactivated FBS | Proliferative phase normal endometrium, 12 days, replicate3 |
| GSM1954917 | Control | treatment: 12d 10% charcoal-stripped heat-inactivated FBS | Proliferative phase normal endometrium, 12 days, replicate4 |
| GSM1954918 | Case | treatment: 12d dibutyryl-cAMP and dienogest | Proliferative phase normal endometrium, 12 days, dibutyryl-cAMP and dienogest, replicate1 |
| GSM1954919 | Case | treatment: 12d dibutyryl-cAMP and dienogest | Proliferative phase normal endometrium, 12 days, dibutyryl-cAMP and dienogest, replicate2 |
| GSM1954920 | Case | treatment: 12d dibutyryl-cAMP and dienogest | Proliferative phase normal endometrium, 12 days, dibutyryl-cAMP and dienogest, replicate3 |
| GSM1954921 | Case | treatment: 12d dibutyryl-cAMP and dienogest | Proliferative phase normal endometrium, 12 days, dibutyryl-cAMP and dienogest, replicate4 |

### Sample Sheet for GSE51981

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM1256653 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256654 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256655 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256656 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256657 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256658 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256659 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256660 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256661 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256662 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256663 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256664 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256665 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256666 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256667 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256668 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256669 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256670 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256671 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256672 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256673 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256674 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256675 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256676 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256677 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256678 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256679 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256680 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256681 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256682 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256683 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256684 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256685 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256686 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256687 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256688 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256689 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256690 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256691 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256692 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256693 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256694 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256695 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256696 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256697 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256698 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256699 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256700 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256701 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256702 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256703 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256704 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256705 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256706 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256707 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256708 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256709 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256710 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256711 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256712 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256713 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256714 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256715 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256716 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256717 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256718 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256719 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256720 | Control | Control | Proliferative Endometrial tissue |
| GSM1256721 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256722 | Control | Control | Proliferative Endometrial tissue |
| GSM1256723 | Control | Control | Proliferative Endometrial tissue |
| GSM1256724 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256725 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256726 | Control | Control | Proliferative Endometrial tissue |
| GSM1256727 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256728 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256729 | Control | Control | Proliferative Endometrial tissue |
| GSM1256730 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256731 | Control | Control | Proliferative Endometrial tissue |
| GSM1256732 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256733 | Control | Control | Proliferative Endometrial tissue |
| GSM1256734 | Control | Control | Proliferative Endometrial tissue |
| GSM1256735 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256736 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256737 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256738 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256739 | Control | Control | Proliferative Endometrial tissue |
| GSM1256740 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256741 | Control | Control | Proliferative Endometrial tissue |
| GSM1256742 | Control | Control | Proliferative Endometrial tissue |
| GSM1256743 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256744 | Control | Control | Proliferative Endometrial tissue |
| GSM1256745 | Control | Control | Proliferative Endometrial tissue |
| GSM1256746 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256747 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256748 | Control | Control | Proliferative Endometrial tissue |
| GSM1256749 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256750 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256751 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256752 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256753 | Control | Control | Proliferative Endometrial tissue |
| GSM1256754 | Control | Control | Proliferative Endometrial tissue |
| GSM1256755 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256756 | Control | Control | Proliferative Endometrial tissue |
| GSM1256757 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256758 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256759 | Control | Control | Proliferative Endometrial tissue |
| GSM1256760 | Control | Control | Proliferative Endometrial tissue |
| GSM1256761 | Control | Control | Proliferative Endometrial tissue |
| GSM1256762 | Control | Control | Proliferative Endometrial tissue |
| GSM1256763 | Control | Control | Proliferative Endometrial tissue |
| GSM1256764 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256765 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256766 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256767 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256768 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256769 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256770 | Control | Control | Proliferative Endometrial tissue |
| GSM1256771 | Control | Control | Proliferative Endometrial tissue |
| GSM1256772 | Control | Control | Proliferative Endometrial tissue |
| GSM1256773 | Case | Endometriosis | Unknown |
| GSM1256774 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256775 | Case | Endometriosis | Unknown |
| GSM1256776 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256777 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256778 | Case | Endometriosis | Early Secretory Endometrial tissue |
| GSM1256779 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256780 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256781 | Case | Endometriosis | Mid-Secretory Endometrial tissue |
| GSM1256782 | Case | Endometriosis | Proliferative Endometrial tissue |
| GSM1256783 | Control | Control | Late Secretory Endometrial tissue |
| GSM1256784 | Control | Control | Proliferative Endometrial tissue |
| GSM1256785 | Control | Control | Proliferative Endometrial tissue |
| GSM1256786 | Control | Control | Proliferative Endometrial tissue |
| GSM1256787 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256788 | Control | Control | Proliferative Endometrial tissue |
| GSM1256789 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256790 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256791 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256792 | Control | Control | Early Secretory Endometrial tissue |
| GSM1256793 | Control | Control | Proliferative Endometrial tissue |
| GSM1256794 | Control | Control | Proliferative Endometrial tissue |
| GSM1256795 | Control | Control | Late Secretory Endometrial tissue |
| GSM1256796 | Control | Control | Proliferative Endometrial tissue |
| GSM1256797 | Control | Control | Mid-Secretory Endometrial tissue |
| GSM1256798 | Control | Control | Proliferative Endometrial tissue |
| GSM1256799 | Control | Control | Proliferative Endometrial tissue |
| GSM1256800 | Control | Control | Proliferative Endometrial tissue |

### Sample Sheet for GSE6364

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM150190 | Case | Endometriosis | uterus - Proliferative endometrium from patient with  endometriosis |
| GSM150191 | Case | Endometriosis | uterus - Proliferative endometrium from patient with  endometriosis |
| GSM150192 | Case | Endometriosis | uterus - Proliferative endometrium from patient with  endometriosis |
| GSM150193 | Case | Endometriosis | uterus - Proliferative endometrium from patient with  endometriosis |
| GSM150194 | Case | Endometriosis | uterus - Proliferative endometrium from patient with  endometriosis |
| GSM150195 | Case | Endometriosis | uterus - Proliferative endometrium from patient with  endometriosis |
| GSM150196 | Control | Control | uterus - Proliferative endometrium from normal patient |
| GSM150197 | Control | Control | uterus - Proliferative endometrium from normal patient |
| GSM150198 | Control | Control | uterus - Proliferative endometrium from normal patient |
| GSM150199 | Control | Control | uterus - Proliferative endometrium from normal patient |
| GSM150201 | Control | Control | uterus - Proliferative endometrium from normal patient |
| GSM150202 | Case | Endometriosis | uterus - Early secretory endometrium from patient with endometriosis |
| GSM150203 | Case | Endometriosis | uterus - Early secretory endometrium from patient with endometriosis |
| GSM150204 | Case | Endometriosis | uterus - Early secretory endometrium from patient with endometriosis |
| GSM150205 | Case | Endometriosis | uterus - Early secretory endometrium from patient with endometriosis |
| GSM150206 | Case | Endometriosis | uterus - Early secretory endometrium from patient with endometriosis |
| GSM150207 | Case | Endometriosis | uterus - Early secretory endometrium from patient with endometriosis |
| GSM150208 | Control | Control | uterus - Early secretory endometrium from normal patient |
| GSM150209 | Control | Control | uterus - Early secretory endometrium from normal patient |
| GSM150210 | Control | Control | uterus - Early secretory endometrium from normal patient |
| GSM150211 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150212 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150213 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150214 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150215 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150216 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150217 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150218 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150219 | Case | Endometriosis | uterus - Mid secretory endometrium from patient with  endometriosis |
| GSM150220 | Control | Control | uterus - Mid secretory endometrium from normal patient |
| GSM150221 | Control | Control | uterus - Mid secretory endometrium from normal patient |
| GSM150222 | Control | Control | uterus - Mid secretory endometrium from normal patient |
| GSM150223 | Control | Control | uterus - Mid secretory endometrium from normal patient |
| GSM150224 | Control | Control | uterus - Mid secretory endometrium from normal patient |
| GSM150225 | Control | Control | uterus - Mid secretory endometrium from normal patient |
| GSM150226 | Control | Control | uterus - Mid secretory endometrium from normal patient |
| GSM150227 | Control | Control | uterus - Mid secretory endometrium from normal patient |

### Sample Sheet for GSE286315

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM8723559 | Control | Control EcSCs | Endometriosis tissue |
| GSM8723556 | Case | HOXC4 KO | Endometriosis tissue |
| GSM8723560 | Control | Control EcSCs | Endometriosis tissue |
| GSM8723557 | Case | HOXC4 KO | Endometriosis tissue |
| GSM8723561 | Control | Control EcSCs | Endometriosis tissue |
| GSM8723558 | Case | HOXC4 KO | Endometriosis tissue |

### Sample Sheet for GSE102293

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM2733413 | Control | phenotype: normal | ovarian granulosa cells |
| GSM2733414 | Control | phenotype: normal | ovarian granulosa cells |
| GSM2733415 | Control | phenotype: normal | ovarian granulosa cells |
| GSM2733416 | Control | phenotype: normal | ovarian granulosa cells |
| GSM2733417 | Case | phenotype: PCOS (Polycystic ovary syndrome) | ovarian granulosa cells |
| GSM2733418 | Case | phenotype: PCOS (Polycystic ovary syndrome) | ovarian granulosa cells |

### Sample Sheet for GSE106724

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM2849384 | Case | PCOS with hyperandrogenism | Ovarian granulosa cells, PCOS with hyperandrogenism |
| GSM2849385 | Case | PCOS with hyperandrogenism | Ovarian granulosa cells, PCOS with hyperandrogenism |
| GSM2849386 | Case | PCOS with hyperandrogenism | Ovarian granulosa cells, PCOS with hyperandrogenism |
| GSM2849388 | Control | normal control | Ovarian granulosa cells, normal control |
| GSM2849391 | Control | normal control | Ovarian granulosa cells, normal control |
| GSM2849392 | Control | normal control | Ovarian granulosa cells, normal control |
| GSM2849393 | Case | PCOS with hyperandrogenism | Ovarian granulosa cells, PCOS with hyperandrogenism |
| GSM2849395 | Control | normal control | Ovarian granulosa cells, normal control |

### Sample Sheet for GSE114419

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM3141461 | Control | healthy control | granulosa cells from normal patients |
| GSM3141462 | Control | healthy control | granulosa cells from normal patients |
| GSM3141463 | Control | healthy control | granulosa cells from normal patients |
| GSM3141464 | Case | Polycystic Ovary Syndrome (PCOS) | granulosa cells from PCOS patients |
| GSM3141465 | Case | Polycystic Ovary Syndrome (PCOS) | granulosa cells from PCOS patients |
| GSM3141466 | Case | Polycystic Ovary Syndrome (PCOS) | granulosa cells from PCOS patients |

### Sample Sheet for GSE137684

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM4084711 | Case | condition: Normoandrogenic PCOS | ovary |
| GSM4084712 | Control | condition: Normal | ovary |
| GSM4084716 | Control | condition: Normal | ovary |
| GSM4084718 | Control | condition: Normal | ovary |
| GSM4084719 | Control | condition: Normal | ovary |
| GSM4084720 | Case | condition: Normoandrogenic PCOS | ovary |
| GSM4084721 | Case | condition: Normoandrogenic PCOS | ovary |
| GSM4084722 | Case | condition: Normoandrogenic PCOS | ovary |

### Sample Sheet for GSE34526

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM850527 | Control | Control | Normal, biological rep 1 |
| GSM850528 | Control | Control | Normal, biological rep 2 |
| GSM850529 | Control | Control | Normal, biological rep 3 |
| GSM850530 | Case | PCOS | PCOS, biological rep 1 |
| GSM850531 | Case | PCOS | PCOS, biological rep 2 |
| GSM850532 | Case | PCOS | PCOS, biological rep 3 |
| GSM850533 | Case | PCOS | PCOS, biological rep 4 |
| GSM850534 | Case | PCOS | PCOS, biological rep 5 |
| GSM850535 | Case | PCOS | PCOS, biological rep 6 |
| GSM850536 | Case | PCOS | PCOS, biological rep 7 |

### Sample Sheet for GSE6798

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM155631 | Control | Control | Vastus lateralis muscle |
| GSM155643 | Control | Control | Vastus lateralis muscle |
| GSM155644 | Control | Control | Vastus lateralis muscle |
| GSM155729 | Control | Control | Vastus lateralis muscle |
| GSM156170 | Control | Control | Vastus lateralis muscle |
| GSM156171 | Control | Control | Vastus lateralis muscle |
| GSM156176 | Control | Control | Vastus lateralis muscle |
| GSM156177 | Control | Control | Vastus lateralis muscle |
| GSM156178 | Control | Control | Vastus lateralis muscle |
| GSM156179 | Control | Control | Vastus lateralis muscle |
| GSM156180 | Control | Control | Vastus lateralis muscle |
| GSM156181 | Control | Control | Vastus lateralis muscle |
| GSM156184 | Control | Control | Vastus lateralis muscle |
| GSM156186 | Case | PCOS | Vastus lateralis muscle |
| GSM156187 | Case | PCOS | Vastus lateralis muscle |
| GSM156510 | Case | PCOS | Vastus lateralis muscle |
| GSM156511 | Case | PCOS | Vastus lateralis muscle |
| GSM156512 | Case | PCOS | Vastus lateralis muscle |
| GSM156749 | Case | PCOS | Vastus lateralis muscle |
| GSM156750 | Case | PCOS | Vastus lateralis muscle |
| GSM156751 | Case | PCOS | Vastus lateralis muscle |
| GSM156752 | Case | PCOS | Vastus lateralis muscle |
| GSM156753 | Case | PCOS | Vastus lateralis muscle |
| GSM156763 | Case | PCOS | Vastus lateralis muscle |
| GSM156946 | Case | PCOS | Vastus lateralis muscle |
| GSM156948 | Case | PCOS | Vastus lateralis muscle |
| GSM156949 | Case | PCOS | Vastus lateralis muscle |
| GSM156950 | Case | PCOS | Vastus lateralis muscle |
| GSM156951 | Case | PCOS | Vastus lateralis muscle |

### Sample Sheet for GSE8157

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM201849 | Control | Control | Vastus lateralis muscle |
| GSM201850 | Control | Control | Vastus lateralis muscle |
| GSM201851 | Control | Control | Vastus lateralis muscle |
| GSM201852 | Control | Control | Vastus lateralis muscle |
| GSM201853 | Control | Control | Vastus lateralis muscle |
| GSM201854 | Control | Control | Vastus lateralis muscle |
| GSM201855 | Control | Control | Vastus lateralis muscle |
| GSM201856 | Control | Control | Vastus lateralis muscle |
| GSM201857 | Control | Control | Vastus lateralis muscle |
| GSM201858 | Control | Control | Vastus lateralis muscle |
| GSM201859 | Control | Control | Vastus lateralis muscle |
| GSM201861 | Control | Control | Vastus lateralis muscle |
| GSM201862 | Control | Control | Vastus lateralis muscle |
| GSM201863 | Case | PCOS | Vastus lateralis muscle |
| GSM201864 | Case | PCOS | Vastus lateralis muscle |
| GSM201865 | Case | PCOS | Vastus lateralis muscle |
| GSM201866 | Case | PCOS | Vastus lateralis muscle |
| GSM201867 | Case | PCOS | Vastus lateralis muscle |
| GSM201868 | Case | PCOS | Vastus lateralis muscle |
| GSM201869 | Case | PCOS | Vastus lateralis muscle |
| GSM201870 | Case | PCOS | Vastus lateralis muscle |
| GSM201871 | Case | PCOS | Vastus lateralis muscle |
| GSM201872 | Case | PCOS | Vastus lateralis muscle |

### Sample Sheet for GSE199225

| Sample ID | Mapped Group | Original Group Value | Tissue |
| :--- | :--- | :--- | :--- |
| GSM5967063 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967064 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967065 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967066 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967067 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967068 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967069 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967070 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967071 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967072 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967073 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967074 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967075 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967076 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967077 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967078 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967079 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967080 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967081 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967082 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967083 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967084 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967085 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967086 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967087 | Control | Control | myotubes obtained from vastus lateralis |
| GSM5967094 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967095 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967096 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967097 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967098 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967099 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967100 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967101 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967102 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967103 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967104 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967105 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967106 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967107 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967108 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967109 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967110 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967111 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967113 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967114 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967115 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967116 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967117 | Case | PCOS | myotubes obtained from vastus lateralis |
| GSM5967118 | Case | PCOS | myotubes obtained from vastus lateralis |

