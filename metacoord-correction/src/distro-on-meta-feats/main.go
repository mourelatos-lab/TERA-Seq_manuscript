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
	"github.com/jmoiron/sqlx"
)

// Opts is the struct with the options that the program accepts.
type Opts struct {
	DB       string `arg:"required,help:SQLite3 database"`
	Table    string `arg:"required,help:database table name"`
	Where    string `arg:"help:SQL filter injected in WHERE clause"`
	Bed      string `arg:"required,help:BED file with features (at least 6 columns)"`
	Pos      string `arg:"required,help:reference point for reads; one of 5p or 3p"`
	Bins     int    `arg:"required,help:number of bins for each feature"`
	QArea    int    `arg:"--q-area,help:area (nts) around feats to extend query; use a big number."`
	Collapse bool   `arg:"help:collapse reads that have the same pos."`
	ByFeat   bool   `arg:"--by-feat,help:print output by feature."`
	MinLen   int    `arg:"--min-len,help:minimum feature length allowed."`
}

// Version returns the program version.
func (Opts) Version() string {
	return "distro-on-meta-feats 0.2"
}

// Description returns an extended description of the program.
func (Opts) Description() string {
	return "Prints the number of reads in each bin of each feature that is read by a BED file. Each feature is divided into a given number of bins and the number of reads in each bin are counted."
}

func main() {
	var opts Opts
	opts.QArea = 1000
	p := arg.MustParse(&opts)
	if opts.Pos != "5p" && opts.Pos != "3p" {
		p.Fail("--pos must be either 5p or 3p")
	}

	// Get position extracting function
	extractPos := Head
	if opts.Pos == "3p" {
		extractPos = Tail
	}

	// open BED6 scanner
	bedF, err := os.Open(opts.Bed)
	if err != nil {
		log.Fatal(err)
	}
	bedR, err := bed.NewReader(bedF, 6)
	if err != nil {
		log.Fatal(err)
	}
	bedS := featio.NewScanner(bedR)

	// Connect to the database.
	db, err := sqlx.Connect("sqlite3", opts.DB)
	if err != nil {
		log.Fatal(err)
	}

	// Prepare squirrel select builder.
	sel := squirrel.Select("start", "stop").From(opts.Table).Where(
		"strand = ? AND rname = ? AND start BETWEEN ? AND ? AND stop BETWEEN ? AND ?")
	if opts.Where != "" {
		sel = sel.Where(opts.Where)
	}

	// Prepare query
	query, _, err := sel.ToSql()
	if err != nil {
		log.Fatal(err)
	}

	// Prepare range query statement.
	stmt, err := db.Preparex(query)
	if err != nil {
		log.Fatal(err)
	}

	// Loop on the BED entries and measure relative positions.
	switch opts.ByFeat {
	case true:
		fmt.Printf("feat\tbin\tcount\tbinLen\n")
	case false:
		fmt.Printf("bin\tcount\n")
	}
	hist := make([]int, opts.Bins)
	for bedS.Next() {
		b := bedS.Feat()
		toCollapse := make(map[int]bool)

		o, ok := b.(feat.Orienter)
		if !ok {
			log.Fatal("bed feature does not have defined orientation")
		}

		binLen := float64(b.Len()) / float64(opts.Bins)
		ori := o.Orientation()
		start := b.Start() - opts.QArea
		stop := b.End() - 1 + opts.QArea
		rname := b.Location().Name()

		if b.Len() < opts.MinLen {
			if opts.ByFeat {
				for i := 0; i < opts.Bins; i++ {
					fmt.Printf("%s\t%d\t%s\t%s\n", rname, i, "NA", "NA")
					hist[i] = 0
				}
			}
			continue
		}

		// Run the query.
		rows, err := stmt.Queryx(ori, rname, start, stop, start, stop)
		if err != nil {
			log.Fatal(err)
		}

		// Loop on the reads.
		r := &Feature{}
		for rows.Next() {
			if err = rows.StructScan(r); err != nil {
				log.Fatal(err)
			}

			pos := extractPos(r, ori) - b.Start()
			if pos < 0 || pos >= b.Len() {
				continue
			}

			if ori == feat.Reverse {
				pos = b.Len() - pos - 1
			}

			if opts.Collapse {
				if toCollapse[pos] {
					continue
				}
				toCollapse[pos] = true
			}

			bin := int(float64(pos) / binLen)
			if bin < 0 || bin >= opts.Bins {
				log.Fatal("bin out of bounds for ", pos)
			}
			hist[bin]++
		}

		// Print counts if --by-feat is set.
		if opts.ByFeat {
			for i := 0; i < opts.Bins; i++ {
				fmt.Printf("%s\t%d\t%d\t%f\n", rname, i, hist[i], binLen)
				hist[i] = 0
			}
		}

	}
	// Print counts if --by-feat is not set.
	if !opts.ByFeat {
		for i := 0; i < opts.Bins; i++ {
			fmt.Printf("%d\t%d\n", i, hist[i])
		}
	}
}

// Feature is part of an htsdb record that wraps Range and the name of the
// reference.
type Feature struct {
	StartPos int              `db:"start"`
	StopPos  int              `db:"stop"`
	Orient   feat.Orientation `db:"strand"`
}

// Name returns an empty string.
func (e *Feature) Name() string { return "" }

// Start returns the start position of Range.
func (e *Feature) Start() int { return e.StartPos }

// End returns the end position of Feature.
func (e *Feature) End() int { return e.StopPos + 1 }

// Len returns the length of Feature.
func (e *Feature) Len() int { return e.End() - e.Start() }

// Description returns an empty string.
func (e *Feature) Description() string { return "" }

// Location returns the location of Feature.
func (e *Feature) Location() feat.Feature { return nil }

// Orientation returns the orientation of OrientedFeature.
func (e *Feature) Orientation() feat.Orientation { return e.Orient }

// Head returns the head coordinate of r depending on orientation.
func Head(r feat.Range, o feat.Orientation) int {
	if o == feat.Forward {
		return r.Start()
	} else if o == feat.Reverse {
		return r.End() - 1
	}
	panic("orientation must be forward or reverse")
}

// Tail returns the tail coordinate of r depending on orientation.
func Tail(r feat.Range, o feat.Orientation) int {
	if o == feat.Forward {
		return r.End() - 1
	} else if o == feat.Reverse {
		return r.Start()
	}
	panic("orientation must be forward or reverse")
}
