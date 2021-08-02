# Samples
To be able to fully re-analyze the data, please the full TERA-Seq dataset from SRA BioProject accession: [PRJNA673166](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA673166), put the downloaded *fastq.gz* files to directories names by `Library name/fastq` (please see the table bellow) in this directory and rename the downloaded *fastq.gz* files to `reads.1.fastq.gz`.

For example, the prepared Cap-Poly(A) *fastq.gz* have to be placed in `hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/fastq/reads.1.fastq.gz` in this directory. The analysis will not work without the samples directory populated exactly this way!

For basecalling and analysis of poly(A) tails you need to download the corresponding fast5 files and place them in `fast5` directory for each of the libraries. Please note the provided fast5 files are already basecalled (`Guppy 3.3.2-1`). Please note the basecalling was run **without** RTA trimming. For more details on basecalling please see [`run_Guppy.sh`](run_Guppy.sh).

After you prepared the input samples/files you can run `run.sh` script to process the sample files.

| Library type | Sample name | fast5 | fastq | SRA | Library name |
|--|--|--|--|--|--|
| 5TERA | Cap-Poly(A) | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413813) | hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1 |
| 5TERA | Cap & 5P-Poly(A) | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413814) | hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1 |
| 5TERA | 5P-Poly(A) | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.polyA.REL5.long.1.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413817) | hsa.dRNASeq.HeLa.polyA.REL5.long.1 |
| 5TERA | 5OH-Poly(A) | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.polyA.REL5OH.long.1.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX10436130) | hsa.dRNASeq.HeLa.polyA.REL5OH.long.1 |
| 5TERA-short | 5P short-Poly(A) | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.polyA.REL5.1.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413824) | hsa.dRNASeq.HeLa.polyA.REL5.1 |
| 5TERA-short | 5 short-PNK-Poly(A) | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.polyA.PNK.REL5.1.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413815) | hsa.dRNASeq.HeLa.polyA.PNK.REL5.1 |
| TERA3 | TERA3_replicate1 | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.total.REL3.1.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413818) | hsa.dRNASeq.HeLa.total.REL3.1 |
| TERA3 | TERA3_replicate2 | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.total.REL3.2.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413819) | hsa.dRNASeq.HeLa.total.REL3.2 |
| TERA3 | TERA3_replicate3 | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.total.REL3.3.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413820) | hsa.dRNASeq.HeLa.total.REL3.3 |
| 5TERA3 | 5TERA3_replicate1 | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.total.REL5.long.REL3.4.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413821) | hsa.dRNASeq.HeLa.total.REL5.long.REL3.4 |
| 5TERA3 | 5TERA3_replicate2 | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.total.REL5.long.REL3.5.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413822) | hsa.dRNASeq.HeLa.total.REL5.long.REL3.5 |
| 5TERA3 | 5TERA3_replicate3 | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.total.REL5.long.REL3.6.tar.gz) | [download]() | [link](https://www.ncbi.nlm.nih.gov/sra/SRX9413823) | hsa.dRNASeq.HeLa.total.REL5.long.REL3.6 |
| Standard dRNA | CTRL-poly(A) | [download](http://mourelatos.med.upenn.edu/teraseq/fast5/hsa.dRNASeq.HeLa.polyA.1.tar.gz) | [download](https://sra-download.ncbi.nlm.nih.gov/traces/sra74/SRZ/014294/SRR14294791/hsa.dRNASeq.HeLa.polyA.1.fastq.gz) | [link](https://www.ncbi.nlm.nih.gov/sra/SRX10652701) | hsa.dRNASeq.HeLa.polyA.1 |

To re-analyze analysis using external data please populate the additional directories. Follow the same instructions as for TERA-Seq using the table bellow.

| Library type | Sample name | DB id | Library name |
|--|--|--|--|
| Mouse SIRV | SIRV_replicate1 | [ERR2680375](https://www.ebi.ac.uk/ena/browser/view/ERR2680375) | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680375.1 |
| Mouse SIRV | SIRV_replicate2 | [ERR2680379](https://www.ebi.ac.uk/ena/browser/view/ERR2680379) | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680379.1 |
| Mouse SIRV | SIRV_replicate3 | [ERR3363657](https://www.ebi.ac.uk/ena/browser/view/ERR3363657) | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363657.1 |
| Mouse SIRV | SIRV_replicate4 | [ERR3363659](https://www.ebi.ac.uk/ena/browser/view/ERR3363659) | mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363659.1 |
| RNASeq | Illumina (RNA-Seq) | [ENCSR000CPR](https://www.encodeproject.org/files/ENCFF000FOM/) | hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1 |
| RiboSeq | Ribo-Seq | [SRR3306589](https://www.ebi.ac.uk/ena/browser/view/SRR3306589) | hsa.RiboSeq.HeLa.async.2 |
| Akron5Seq | Akron5 | [SRR6360508](https://www.ebi.ac.uk/ena/browser/view/SRR6360508) | hsa.Akron5Seq.HeLa.whole.2 |

> Written with [StackEdit](https://stackedit.io/).
