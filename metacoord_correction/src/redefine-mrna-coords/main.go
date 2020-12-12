package main

import (
	"fmt"
	"log"
	"os"

	_ "github.com/mattn/go-sqlite3"

	"github.com/Masterminds/squirrel"
	"github.com/alexflint/go-arg"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/brentp/xopen"
	"github.com/jmoiron/sqlx"
)

// MaxInt is the maximum integer.
const MaxInt = int(^uint(0) >> 1)

// Opts is the struct with the options that the program accepts.
type Opts struct {
	DB     []string `arg:"required,help:SQLite3 database with reads defining mRNA ends"`
	Table  string   `arg:"required,help:database table name"`
	Where  string   `arg:"help:SQL filter injected in WHERE clause"`
	Bed    string   `arg:"required,help:BED file with mRNA features (at least 6 columns)"`
	BedCDS string   `arg:"required,--bed-cds,help:BED file with the CDS coords on the mRNAs"`
	Thres  int      `arg:"required,help:minimum number of supporting reads for each coordinate"`
}

// Version returns the program version.
func (Opts) Version() string {
	return "redefine-mrna-coords 0.2"
}

// Description returns an extended description of the program.
func (Opts) Description() string {
	return "Given alignments on features it attempts to guess the 5' and 3' ends of the features. It will report the 5' ends and 3' ends that are outside the CDS and are supported by the maximum amount of alignments."
}

func main() {
	var opts Opts
	opts.Thres = 1
	arg.MustParse(&opts)

	if opts.Thres < 1 {
		fmt.Println("error: thres cannot be less than 1")
		os.Exit(1)
	}

	// Open BED6 scanner for mRNAs.
	bedScan, err := bedScanner(opts.Bed)
	if err != nil {
		log.Fatal(err)
	}

	// Open BED6 scanner for CDS regions.
	cdsScan, err := bedScanner(opts.BedCDS)
	if err != nil {
		log.Fatal(err)
	}

	// Store CDS regions in a lookup map based on mRNA name.
	cdss := make(map[string]feat.Feature)
	for cdsScan.Next() {
		b := cdsScan.Feat()
		cdss[b.Location().Name()] = b
	}
	if err = cdsScan.Error(); err != nil {
		log.Fatal(err)
	}

	// Connect to the databases.
	dbs := make([]*sqlx.DB, len(opts.DB))
	for i, dbFile := range opts.DB {
		db, err := sqlx.Connect("sqlite3", dbFile)
		if err != nil {
			log.Fatal(err)
		}
		dbs[i] = db
	}

	// Prepare squirrel select builder.
	sel := squirrel.Select("qname", "start", "stop").From(opts.Table).Where(
		"strand = ? AND rname = ?")
	if opts.Where != "" {
		sel = sel.Where(opts.Where)
	}

	// Prepare query
	query, _, err := sel.ToSql()
	if err != nil {
		log.Fatal(err)
	}

	// Prepare range query statement.
	stmts := make([]*sqlx.Stmt, len(dbs))
	for i, db := range dbs {
		stmt, err := db.Preparex(query)
		if err != nil {
			log.Fatal(err)
		}
		stmts[i] = stmt
	}

	// Loop on the BED mRNA entries.
	for bedScan.Next() {
		b := bedScan.Feat()

		// Features must be on positive strand.
		o, ok := b.(feat.Orienter)
		if !ok {
			log.Fatal("orientation not defined for %v", b)
		}
		if o.Orientation() != feat.Forward {
			log.Fatalf("error: non-forward orientation for %v", b)
		}

		// Get the CDS region.
		cds, ok := cdss[b.Location().Name()]
		if !ok {
			log.Fatalf("error: cds region not found for %v", b)
		}

		// These maps will hold the count for all the start and end positions
		// of the reads on the mRNA.
		starts := make(map[int]int)
		ends := make(map[int]int)

		for _, stmt := range stmts {
			// Run query to retrieve reads within feature.
			rows, err := stmt.Queryx(feat.Forward, b.Location().Name())
			if err != nil {
				log.Fatal(err)
			}

			// Loop on reads and store their start and end in starts and ends.
			// Only starts and ends outside the CDS bounds are stored.
			r := &Feature{}
			for rows.Next() {
				if err = rows.StructScan(r); err != nil {
					log.Fatalf("error: read scanning from db failed: %v", err)
				}

				// The new mRNA coordinates need to be outside the CDS bounds.
				if r.StartPos < cds.Start() {
					starts[r.StartPos]++
				}
				if r.StopPos+1 > cds.End() {
					ends[r.StopPos+1]++
				}
			}
			if err := rows.Err(); err != nil {
				log.Fatalf("error: db iteration failed: %v", err)
			}
		}

		// Define new coords and ensure enough supporting observations.
		newStart, supportCnt := findMinKeyMaxVal(starts)
		if supportCnt < opts.Thres {
			continue
		}
		newEnd, supportCnt := findMaxKeyMaxVal(ends)
		if supportCnt < opts.Thres {
			continue
		}

		// Print BED entry with the new coordinates.
		fmt.Printf("%s\t%d\t%d\t%s\t%d\t%s\t#%d,%d\n",
			b.Location().Name(),
			newStart,
			newEnd,
			"mrna",
			newEnd-newStart,
			"+",
			starts[newStart],
			ends[newEnd],
		)
	}
	if err = bedScan.Error(); err != nil {
		log.Fatal(err)
	}
}

// Feature is part of an htsdb record that wraps Range and the name of the
// reference.
type Feature struct {
	Qname    string           `db:"qname"`
	StartPos int              `db:"start"`
	StopPos  int              `db:"stop"`
	Orient   feat.Orientation `db:"strand"`
}

// findMinKeyMaxVal returns the key, value pair of p with the highest value.
// If two pairs have equal values it returns the pair with the smallest key.
func findMinKeyMaxVal(p map[int]int) (int, int) {
	maxIdx := MaxInt
	maxVal := 0
	for i, v := range p {
		if (v > maxVal) || (v == maxVal && i < maxIdx) {
			maxVal = v
			maxIdx = i
		}
	}

	return maxIdx, maxVal
}

// findMaxKeyMaxVal returns the key, value pair of p with the highest value.
// If two pairs have equal values it returns the pair with the larger key.
func findMaxKeyMaxVal(p map[int]int) (int, int) {
	maxIdx := 0
	maxVal := 0
	for i, v := range p {
		if (v > maxVal) || (v == maxVal && i > maxIdx) {
			maxVal = v
			maxIdx = i
		}
	}

	return maxIdx, maxVal
}

// findMaxVal returns the maximum value in p.
func findMaxVal(p map[int]int) int {
	max := 0
	for _, v := range p {
		if v > max {
			max = v
		}
	}
	return max
}

func bedScanner(f string) (*featio.Scanner, error) {
	bedF, err := xopen.Ropen(f)
	if err != nil {
		return nil, err
	}
	bedR, err := bed.NewReader(bedF, 6)
	if err != nil {
		return nil, err
	}
	return featio.NewScanner(bedR), nil
}
