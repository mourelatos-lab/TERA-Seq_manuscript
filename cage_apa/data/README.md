
If you have Salmon expression tables for the libraries (quant.sf) you can extract transcript names with non-zero (>0) expression:


```
libraries=(
	"hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
)

SAMPLE_DIR="../../samples"

for i in "${libraries[@]}"; do
        cat $SAMPLE_DIR/$i/counts/quant.sf | \
        tail -n+2 | \
        awk -F'\t' '$4>0 {print $1}'
done | sort | uniq > HeLa_transcripts.txt
```
