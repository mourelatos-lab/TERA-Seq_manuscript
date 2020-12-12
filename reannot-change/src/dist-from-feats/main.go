package main

import (
	"fmt"
	"log"
	"os"

	_ "github.com/mattn/go-sqlite3"

	"github.com/Masterminds/squirrel"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/jmoiron/sqlx"
	"gopkg.in/alecthomas/kingpin.v2"
)

const prog = "dist-from-feats"
const version = "0.2"
const descr = `Prints a table with the aggregate number of reads at each position around the 5' or 3' end of the entries in a BED file.`

var (
	app = kingpin.New(prog, descr)

	dbFile = app.Flag("db", "SQLite file.").
		PlaceHolder("<file>").Required().String()
	table = app.Flag("table", "Database table name.").
		Default("sample").String()
	where = app.Flag("where", "SQL filter injected in WHERE clause.").
		PlaceHolder("<SQL>").String()
	bed6 = app.Flag("bed", "BED file with features.").
		PlaceHolder("<file>").Required().String()
	whichBedPos = app.Flag("bed-pos", "Reference point for BED entries.").
			Required().PlaceHolder("<5p|3p>").Enum("5p", "3p")
	whichPos = app.Flag("pos", "Reference point for database entries.").
			Required().PlaceHolder("<5p|3p>").Enum("5p", "3p")
	spanUp = app.Flag("span-up", "Max upstream position around bed-pos.").
		Default("100").PlaceHolder("<int>").Int()
	spanDown = app.Flag("span-down", "Max downstream position around bed-pos.").
			Default("100").PlaceHolder("<int>").Int()
	qArea = app.Flag("q-area", "Area around bed-pos to do a range query; use a big number.").
		Default("1000").PlaceHolder("<int>").Int()
	anti = app.Flag("anti", "Compare reads on opposite instead of same orientation.").
		Bool()
	collapse = app.Flag("collapse", "Collapse the reads that have the same pos.").
			Bool()
)

func main() {
	app.HelpFlag.Short('h')
	app.Version(version)
	_, err := app.Parse(os.Args[1:])
	if err != nil {
		kingpin.Fatalf("%s", err)
	}

	// Get position extracting function
	extractBedPos := Head
	if *whichBedPos == "3p" {
		extractBedPos = Tail
	}
	extractPos := Head
	if *whichPos == "3p" {
		extractPos = Tail
	}

	// open BED6 scanner
	bedS, err := bed6Scanner(*bed6)
	if err != nil {
		panic(err)
	}

	// Connect to the database.
	db, err := sqlx.Connect("sqlite3", *dbFile)
	if err != nil {
		log.Fatal(err)
	}

	// Prepare squirrel select builder.
	sel := squirrel.Select("start", "stop").From(*table).Where(
		"strand = ? AND rname = ? AND start BETWEEN ? AND ? AND stop BETWEEN ? AND ?")

	// Apply extra where clause if provided.
	if *where != "" {
		sel = sel.Where(*where)
	}

	// Prepare range query statement.
	stmt, err := prepareStmt(sel, db)
	if err != nil {
		log.Fatal(err)
	}

	// Loop on the BED entries and measure relative positions.
	fmt.Printf("ref\tfrom\tpos\tcount\n")
	for {
		if bedS.Next() == false {
			break
		}
		b := bedS.Feat()
		hist := make(map[int]int)
		toCollapse := make(map[int]bool)

		o, ok := b.(feat.Orienter)
		if !ok {
			log.Fatal("bed feature does not have defined orientation")
		}
		bori := o.Orientation()

		start := b.Start() - *qArea
		stop := b.End() - 1 + *qArea
		rname := b.Location().Name()
		ori := bori
		if *anti == true {
			ori = -1 * bori
		}

		var rr []Feature
		if err := stmt.Select(&rr, ori, rname, start, stop, start, stop); err != nil {
			log.Fatal(err)
		}

		bPos := extractBedPos(b, bori)
		for _, r := range rr {
			pos := extractPos(&r, ori)
			if *collapse {
				if toCollapse[pos] {
					continue
				}
				toCollapse[pos] = true
			}
			dist := (pos - bPos) * int(bori)
			if dist >= -*spanUp && dist <= *spanDown {
				hist[dist]++
			}
		}

		for i := -*spanUp; i <= *spanDown; i++ {
			fmt.Printf("%s\t%d\t%d\t%d\n", rname, bPos, i, hist[i])
		}
	}
}

// prepareStmt returns an *sqlx.Stmt prepared from the squirrel Builder b.
func prepareStmt(b squirrel.SelectBuilder, db *sqlx.DB) (*sqlx.Stmt, error) {
	q, _, err := b.ToSql()
	if err != nil {
		return nil, err
	}
	stmt, err := db.Preparex(q)
	if err != nil {
		return nil, err
	}
	return stmt, nil
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
func (e *Feature) Location() feat.Feature {
	return nil
}

// Orientation returns the orientation of OrientedFeature.
func (e *Feature) Orientation() feat.Orientation {
	return e.Orient
}

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

func bed6Scanner(f string) (*featio.Scanner, error) {
	ioR, err := os.Open(*bed6)
	if err != nil {
		return nil, err
	}
	bedR, err := bed.NewReader(ioR, 6)
	if err != nil {
		return nil, err
	}
	return featio.NewScanner(bedR), nil
}
