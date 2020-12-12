package main

import (
	"fmt"
	"log"
	"os"
	"sync"

	_ "github.com/mattn/go-sqlite3"

	"github.com/Masterminds/squirrel"
	"github.com/alexflint/go-arg"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/jmoiron/sqlx"
	"github.com/mnsmar/htsdb"
)

const maxConc = 12

// Opts is the struct with the options that the program accepts.
type Opts struct {
	DB     string `arg:"required,help:SQLite3 database"`
	Table  string `arg:"required,help:table name for db"`
	Where  string `arg:"help:SQL filter injected in WHERE clause of db"`
	BED    string `arg:"required,help:BED file with feats (at least 6 columns)"`
	CDS    string `arg:"required,help:BED file with coding regions in feats"`
	RefLim int    `arg:"help:Minimum size of BED entry"`
}

// Version returns the program version.
func (Opts) Version() string {
	return "count-reads-by-decay-level 0.1"
}

// Description returns an extended description of the program.
func (Opts) Description() string {
	return "Count the number of reads on each transcript. Count reads stratifying them to intact CDS and total."
}

func main() {
	var err error
	var db *sqlx.DB

	// reads command line arguments.
	var opts Opts
	arg.MustParse(&opts)

	// open database connection.
	if db, err = sqlx.Connect("sqlite3", opts.DB); err != nil {
		log.Fatal(err)
	}

	// create a list of decorators for the select statement. Only the basic
	// decorators are appended here. The actual select statement is built
	// individually by each worker.
	decors1 := []BuilderDecorator{Table(opts.Table), Where(opts.Where)}

	// open BED6 scanner for mRNAs.
	bedS, err := makeBEDScanner(opts.BED)
	if err != nil {
		log.Fatal(err)
	}

	// open BED6 scanner for CDS and create a lookup map.
	cdsScan, err := makeBEDScanner(opts.CDS)
	if err != nil {
		log.Fatal(err)
	}
	cdss := make(map[string]feat.Feature)
	for cdsScan.Next() {
		if cdsScan.Error() != nil {
			log.Fatal("error reading cds bed: ", bedS.Error())
		}
		b := cdsScan.Feat()
		cdss[b.Location().Name()] = b
	}
	if err = cdsScan.Error(); err != nil {
		log.Fatal("error reading cds bed: ", err)
	}

	// goroutine that sends each reference as a job to jobs.
	jobs := make(chan job)
	go func() {
		for bedS.Next() {
			if bedS.Error() != nil {
				log.Fatal("error reading bed: ", bedS.Error())
			}
			b := bedS.Feat()

			if b.Len() < opts.RefLim {
				continue
			}

			jobs <- job{
				opts:    opts,
				ref:     b,
				db:      db,
				decors1: decors1,
				cds:     cdss[b.Location().Name()],
			}
		}
		if err = bedS.Error(); err != nil {
			log.Fatal("error reading bed: ", err)
		}
		close(jobs)
	}()

	// start workers that consume jobs and send results to results.
	results := make(chan result)
	var wg sync.WaitGroup
	wg.Add(maxConc)
	for w := 1; w <= maxConc; w++ {
		go func() {
			worker(w, jobs, results)
			wg.Done()
		}()
	}

	// goroutine that checks when all workers are done and closes results.
	go func() {
		wg.Wait()
		close(results)
	}()

	// print output
	fmt.Printf("ref\ttotal\tfull\tfullcds\tlen\tcds_len\n")
	for res := range results {
		if res.hasCDS {
			fmt.Printf("%s\t%d\t%d\t%d\t%d\t%d\n",
				res.job.ref.Location().Name(),
				res.total, res.fullLen, res.fullCDS,
				res.job.ref.Len(), res.job.cds.Len())
		} else {
			fmt.Printf("%s\t%d\t%d\t%s\t%d\t%s\n",
				res.job.ref.Location().Name(),
				res.total, res.fullLen, "NA",
				res.job.ref.Len(), "NA")
		}

	}
}

func worker(id int, jobs <-chan job, results chan<- result) {
	for job := range jobs {
		var err error
		ref := job.ref
		cds := job.cds

		// assemble sqlx select builders
		rangeDec := Where("rname = ?")
		readsB1 := DecorateBuilder(htsdb.RangeBuilder, append(job.decors1, rangeDec)...)

		// prepare statements.
		var readsStmt1 *sqlx.Stmt
		if readsStmt1, err = prepareStmt(readsB1, job.db); err != nil {
			log.Fatal(err)
		}

		var r htsdb.Range
		res := result{
			job:    job,
			hasCDS: cds != nil,
		}

		rows1, err := readsStmt1.Queryx(ref.Location().Name())
		if err != nil {
			log.Fatal(err)
		}
		for rows1.Next() {
			if err = rows1.StructScan(&r); err != nil {
				log.Fatal(err)
			}
			res.total++

			if cds != nil {
				if r.Start() <= cds.Start() && r.End() >= cds.End() {
					res.fullCDS++
				}
			}

			margin := int(float64(ref.Len()) * 0.05) // 50nts for 1000nt mRNA
			if r.Start() <= ref.Start()+margin && r.End() >= ref.End()-margin {
				res.fullLen++
			}
		}

		// enqueue in results channel
		results <- res
	}
}

type job struct {
	opts     Opts
	ref, cds feat.Feature
	decors1  []BuilderDecorator
	db       *sqlx.DB
}

type result struct {
	total   int
	fullCDS int
	fullLen int
	job     job
	hasCDS  bool
}

// A BuilderDecorator wraps a squirrel.SelectBuilder with extra behaviour.
type BuilderDecorator func(squirrel.SelectBuilder) squirrel.SelectBuilder

// Table returns a BuilderDecorator that extends a squirrel.SelectBuilder with
// the table property.
func Table(table string) BuilderDecorator {
	return func(b squirrel.SelectBuilder) squirrel.SelectBuilder {
		return b.From(table)
	}
}

// Where returns a BuilderDecorator that extends a squirrel.SelectBuilder with
// a where clause. Returns the builder itself if the where clause is the empty
// string.
func Where(clause string) BuilderDecorator {
	return func(b squirrel.SelectBuilder) squirrel.SelectBuilder {
		if clause != "" {
			return b.Where(clause)
		}
		return b
	}
}

//DecorateBuilder decorates a squirrel.SelectBuilder with all the given
//BuilderDecorators, in order.
func DecorateBuilder(b squirrel.SelectBuilder, ds ...BuilderDecorator) squirrel.SelectBuilder {
	decorated := b
	for _, decorate := range ds {
		decorated = decorate(decorated)
	}
	return decorated
}

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

func makeBEDScanner(filepath string) (*featio.Scanner, error) {
	// open BED6 scanner for mRNAs.
	bedF, err := os.Open(filepath)
	if err != nil {
		return nil, err
	}
	bedR, err := bed.NewReader(bedF, 6)
	if err != nil {
		return nil, err
	}
	return featio.NewScanner(bedR), nil
}
