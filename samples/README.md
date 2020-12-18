# Samples
To be able to fully re-analyze the data, please the full TERA-Seq dataset from SRA BioProject accession: [PRJNA673166](https://www.ncbi.nlm.nih.gov/bioproject/673166), put the downloaded *fastq.gz* files to directories names by `Library name/fastq` (please see the table bellow) in this directory and rename the downloaded *fastq.gz* files to `reads.1.fastq.gz`. 

For example, the prepared Cap-Poly(A) *fastq.gz* have to be placed in `hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/fastq/reads.1.fastq.gz` in this directory. The analysis will not work without the samples directory populated exactly this way!

After you prepared the input samples/files you can run `run.sh` script to process the sample files. 

Please contact [Jan Oppelt](mailto:jan.oppelt@pennmedicine.upenn.edu) if you need to get the corresponding *fast5* files.

| Library type | Sample name| DB id | Library name |
|--|--|--|--|
| 5TERA | Cap-Poly(A) | [TBU]() | hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1 |
| 5TERA | Cap & 5P-Poly(A) | [TBU]() | hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1 |
| 5TERA | 5P-Poly(A) | [TBU]() | hsa.dRNASeq.HeLa.polyA.REL5.long.1 |
| 5TERA-short | 5P short-Poly(A) | [TBU]() | hsa.dRNASeq.HeLa.polyA.REL5.1 |
| 5TERA-short | 5 short-PNK-Poly(A) | [TBU]() | hsa.dRNASeq.HeLa.polyA.PNK.REL5.1 |
| TERA3 | TERA3_replicate1 | [TBU]() | hsa.dRNASeq.HeLa.total.REL3.1 |
| TERA3 | TERA3_replicate2 | [TBU]() | hsa.dRNASeq.HeLa.total.REL3.2 |
| TERA3 | TERA3_replicate3 | [TBU]() | hsa.dRNASeq.HeLa.total.REL3.3 |
| 5TERA3 | 5TERA3 | [TBU]() | hsa.dRNASeq.HeLa.total.REL5.long.REL3.4 |
| 5TERA3 | 5TERA3 | [TBU]() | hsa.dRNASeq.HeLa.total.REL5.long.REL3.5 |
| 5TERA3 | 5TERA3 | [TBU]() | hsa.dRNASeq.HeLa.total.REL5.long.REL3.6 |
| 5TERA | 5TERA SIRV | [TBU]() | hsa.dRNASeq.SIRV.polyA.REL5.long.2 |

Note: The direct link to the TERA-Seq datasets will be updated once the manuscript has been published. Until then, the datasets are available only for reviewers.

To re-analyze analysis using external data please populate the additional directories. Follow the same instructions as for TERA-Seq using the table bellow.

| Library type | Sample name | DB id | Library name |
|--|--|--|--|
| Mouse SIRV | SIRV | [ERR2680375](https://www.ebi.ac.uk/ena/browser/view/ERR2680375) | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680375.1 |
| Mouse SIRV | SIRV | [ERR2680379](https://www.ebi.ac.uk/ena/browser/view/ERR2680379)  | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680379.1 |
| Mouse SIRV | SIRV | [ERR3363657](https://www.ebi.ac.uk/ena/browser/view/ERR3363657)  | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363657.1 |
| Mouse SIRV | SIRV | [ERR3363659](https://www.ebi.ac.uk/ena/browser/view/ERR3363659)  | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363659.1 |
| RNASeq | Illumina (RNA-Seq) | [ENCSR000CPR](https://www.encodeproject.org/files/ENCFF000FOM/) | hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1 |
| RiboSeq | Ribo-Seq | [SRR3306589](https://www.ebi.ac.uk/ena/browser/view/SRR3306589) | hsa.RiboSeq.HeLa.async.2 |
| Akron5Seq | Akron5 | [SRR6360508](https://www.ebi.ac.uk/ena/browser/view/SRR6360508) | hsa.Akron5Seq.HeLa.whole.2 |

> Written with [StackEdit](https://stackedit.io/).
