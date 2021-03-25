package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"regexp"

	_ "github.com/mattn/go-sqlite3"

	"github.com/Masterminds/squirrel"
	"github.com/alexflint/go-arg"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jmoiron/sqlx"
	"github.com/mnsmar/wig/rand"
)

// Opts is the struct with the options that the program accepts.
type Opts struct {
	DB         string `arg:"required,help:SQLite3 database"`
	Table      string `arg:"required,help:database table name"`
	Where      string `arg:"help:SQL filter injected in WHERE clause"`
	JSON       string `arg:"required,help:JSON file with values for a variable along transcripts"`
	Pos        string `arg:"required,help:reference point for reads; one of 5p or 3p"`
	SpanUp     int    `arg:"--span-up,required,help:max upstream distance from pos"`
	SpanDown   int    `arg:"--span-down,required,help:max downstream distance from pos"`
	Collapse   bool   `arg:"help:collapse reads that have the same pos"`
	BED        string `arg:"required,help:BED file with feats (at least 6 columns)"`
	LenLim     int    `arg:"--len-lim,help:minimum feature length"`
	QArea      int    `arg:"--q-area,help:area (nts) to extend query around feats; use a big number Default: 1000"`
	Random     bool   `arg:"help:randomise pos"`
	RandomKeep bool   `arg:"--random-keep,help:enable randomization preserving nt content"`
	FASTA      string `arg:"required,help:FASTA file with reference sequences; enables randomization preserving nt content"`
	WinUp      int    `arg:"--win-up,help:max upstream distance from pos to preserve nt content"`
	WinDown    int    `arg:"--win-down,help:max downstream distance from pos to preserve nt content"`
	NoGquad    bool   `arg:"--no-gquad,help:keep records without downstream g-gquad"`
}

// Version returns the program version.
func (Opts) Version() string {
	return "htsdb-var-distro-with-gquad-in-feats 0.1"
}

// Description returns an extended description of the program.
func (Opts) Description() string {
	return "Prints the distribution of a variable around the 5'/3' end of reads located with features defined in a BED file. The variable values are defined on each position of the reference sequence and are read from a JSON file. Read positions with no sufficient upstream and downstream span space to the feature ends are skipped. Only positions where a G-Quad is found within 100 nts downstream are retained"
}

type values []float64

func main() {
	var opts Opts
	opts.QArea = 1000
	p := arg.MustParse(&opts)
	if opts.Pos != "5p" && opts.Pos != "3p" {
		p.Fail("--pos must be either 5p or 3p")
	}
	if opts.SpanUp < 0 && opts.SpanDown < 0 {
		p.Fail("--span-up and --span-down must be positive integers")
	}

	// Define the G-Quadruplex prediction pattern.
	validG4 := regexp.MustCompile(`G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}`)
	margin := 100 // min sequence length upstream and downstream

	// get position extracting function
	extractPos := Head
	if opts.Pos == "3p" {
		extractPos = Tail
	}

	// Open JSON and store data in dat.
	var dat map[string]values
	raw, err := ioutil.ReadFile(opts.JSON)
	if err != nil {
		panic(err)
	}
	if err := json.Unmarshal(raw, &dat); err != nil {
		panic(err)
	}

	// Read FASTA records in map.
	seqs := make(map[string][]byte)
	if opts.FASTA != "" {
		fastaF, err := os.Open(opts.FASTA)
		if err != nil {
			panic(err)
		}
		fastaR := fasta.NewReader(fastaF, linear.NewSeq("", nil, alphabet.DNA))
		fastaS := seqio.NewScanner(fastaR)
		for fastaS.Next() {
			s := fastaS.Seq().(*linear.Seq)
			seqs[s.Name()] = alphabet.LettersToBytes(s.Seq)
		}
	}

	// open BED6 scanner
	bedF, err := os.Open(opts.BED)
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
	sel := squirrel.Select("strand", "rname", "start", "stop").From(opts.Table).Where(
		"strand = ? AND rname = ? AND start BETWEEN ? AND ? AND stop BETWEEN ? AND ?")

	// Apply extra where clause if provided.
	if opts.Where != "" {
		sel = sel.Where(opts.Where)
	}

	// Prepare SQL query.
	query, _, err := sel.ToSql()
	if err != nil {
		log.Fatal(err)
	}

	// Prepare statement.
	stmt, err := db.Preparex(query)
	if err != nil {
		log.Fatal(err)
	}

	// Loop on the features.
	counts := make(map[int]float64)
	recordCnt := 0
	for bedS.Next() {
		b := bedS.Feat()
		if opts.LenLim != 0 && b.Len() < opts.LenLim {
			continue
		}
		if b.Len() <= opts.SpanUp+opts.SpanDown {
			continue
		}

		o, ok := b.(feat.Orienter)
		if !ok {
			log.Fatal("bed feature does not have defined orientation")
		}

		qOri := o.Orientation()
		if qOri == feat.Reverse {
			log.Fatal("features on reverse orientation are not supported")
		}

		qRname := b.Location().Name()
		vals, ok := dat[qRname]
		if !ok {
			continue
		}

		seq, ok := seqs[qRname]
		if !ok {
			continue
		}

		toCollapse := make(map[int]bool)
		wig := make([]int, b.End())
		qStart := b.Start() - opts.QArea
		qStop := b.End() - 1 + opts.QArea
		lowerPosBound := b.Start() + opts.SpanUp
		upperPosBound := b.End() - opts.SpanDown

		rows, err := stmt.Queryx(qOri, qRname, qStart, qStop, qStart, qStop)
		if err != nil {
			log.Fatal(err)
		}

		r := &Feature{}
		for rows.Next() {
			if err = rows.StructScan(r); err != nil {
				log.Fatal(err)
			}
			pos := extractPos(r, qOri)

			// leave enough space for spans.
			if pos < lowerPosBound || pos >= upperPosBound {
				continue
			}

			if opts.Collapse {
				if toCollapse[pos] {
					continue
				}
				toCollapse[pos] = true
			}

			wig[pos]++
		}

		if opts.Random {
			rand.Perm(wig[lowerPosBound:upperPosBound], nil)
		}
		if opts.RandomKeep {
			for orf := 0; orf < 3; orf++ {
				rand.PermKeepByte(
					wig[lowerPosBound:upperPosBound],
					seq[lowerPosBound:upperPosBound],
					-opts.WinUp, opts.WinDown,
					func(i int) bool {
						if (i+opts.SpanUp)%3 == orf {
							return true
						}
						return false
					})
			}
		}

		for pos := lowerPosBound; pos < upperPosBound; pos++ {
			if wig[pos] < 1 {
				continue
			}

			subseq := seq[pos : pos+margin]
			g4found := validG4.Match(subseq)
			if g4found && opts.NoGquad {
				continue
			} else if !g4found && !opts.NoGquad {
				continue
			}

			for i := -opts.SpanUp; i <= opts.SpanDown; i++ {
				counts[i] += float64(wig[pos]) * vals[pos+i]
			}
			recordCnt += wig[pos]
		}
	}
	if err = bedS.Error(); err != nil {
		log.Fatal(err)
	}

	// Print the counts.
	fmt.Print("pos\tval\trecords\n")
	for pos := -opts.SpanUp; pos <= opts.SpanDown; pos++ {
		fmt.Printf("%d\t%f\t%d\n", pos, counts[pos], recordCnt)
	}
}

// Feature is part of an htsdb record that wraps Range and the name of the
// reference.
type Feature struct {
	Rname    string           `db:"rname"`
	Orient   feat.Orientation `db:"strand"`
	StartPos int              `db:"start"`
	StopPos  int              `db:"stop"`
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
