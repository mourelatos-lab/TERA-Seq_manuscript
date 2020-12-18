# Meta-coordinates, re-annotation and heatmap

## Corrected coordinates
The corrected coordinates can be found in `results/coords-corrected/new-mrna-coords.5p.unambig.tab` for 5' end coordinates and `results/coords-corrected/new-mrna-coords.3p.unambig.tab` for 3' end coordinates. The orignal transcript annotation with **included corrections** can be found in `results/coords-corrected/genic_elements.mrna.unambig.tab`.

## Distribution along meta-coordinates
### All transcripts together
The **meta-coordinates distribution** plots along meta-coordinates with the **default annotation** for 5TERA libraries can be found in `results/library-name/uncorrected/coords-on-meta-mrnas.unambig.rel5.pdf` for reads with 5TERA adapter, in `results/library-name/uncorrected/coords-on-meta-mrnas.unambig.norel5.pdf` for reads without 5TERA adapter, or in `results/library-name/uncorrected/coords-on-meta-mrnas.unambig.norel5_vs_rel5.pdf` for comparison of both w/o 5TERA adapter reads. *library-name* is a placeholder for library name.

The **meta-coordinates distribution** plots along meta-coordinates with the **corrected annotation** can be found in `results/library-name/corrected/coords-on-corrected-meta-mrnas.unambig.rel5.pdf` for reads with 5TERA adapter, in `results/library-name/corrected/coords-on-corrected-meta-mrnas.unambig.norel5.pdf` for reads without 5TERA adapter, or in `results/library-name/corrected/coords-on-corrected-meta-mrnas.unambig.norel5_vs_rel5.pdf` for comparison of both w/o 5TERA adapter reads.

Page 4 contains distribution of reads averaged to transcripts as used in the manuscript. Page 1 contains distribution of reads without transcript averaging (each read separately).

The same naming scheme applies for TERA3 libraries.

### Transcripts categorized by length
Plots, where transcripts are categorized/divided by length can be found in the the same files as in the plots above with **-inclLen-sameQuart.pdf** suffix added to the original file name.

Page 4 contains distribution of reads averaged to transcripts as used in the manuscript. Page 1 contains distribution of reads without transcript averaging (each read separately).

### 5TERA3 distribution w/o poly(A) tail
**Note:** This analysis can be done only if you ran Nanopolish poly(A) estimates during the sample preparation (`samples` directory).

Visualization of 5TERA3 read distribution of reads **with** and **without poly(A) tail** can be found at `results/hsa.dRNASeq.HeLa.total.REL5.long.REL3.X/corrected/coords-on-corrected-meta-mrnas.unambig.rel5.coding.w_polya.pdf` for reads **with** poly(A) tail; and `results/hsa.dRNASeq.HeLa.total.REL5.long.REL3.X/corrected/coords-on-corrected-meta-mrnas.unambig.rel5.coding.wo_polya.pdf` for reads **without** poly(A) tail.

Page 4 contains distribution of reads averaged to transcripts as used in the manuscript. Page 1 contains distribution of reads without transcript averaging (each read separately).

## Heatmap
**Heatmap distribution**  with the **default annotation** can be found in `results/library-name/uncorrected/coords-on-meta-mrnas-heatmap.unambig.rel5.pdf` for reads with 5TERA adapter, or in `results/library-name/uncorrected/coords-on-meta-mrnas-heatmap.unambig.norel5.pdf` for reads without 5TERA adapter.

**Heatmap distribution**  with the **corrected annotation** can be found in `results/library-name/corrected/coords-on-corrected-meta-mrnas-heatmap.unambig.rel5.pdf` for reads with 5TERA adapter, or in `results/library-name/corrected/coords-on-corrected-meta-mrnas-heatmap.unambig.norel5.pdf` for reads without 5TERA adapter.

The same naming scheme applies to TERA3 libraries.

## Note
Alternative corrections and plots using *primary* alignments instead *unambigous* is available as well (look for *primary* instead of *unambig* in the file name).
