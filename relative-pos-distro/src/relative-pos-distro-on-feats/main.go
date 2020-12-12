package main

import (
	"bytes"
	"fmt"
	"log"
	"os"
	"sync"

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
	"github.com/mnsmar/htsdb"
	"github.com/mnsmar/wig/rand"
)

const maxConc = 12

// Opts is the struct with the options that the program accepts.
type Opts struct {
	DB1        string `arg:"required,help:SQLite3 database 1"`
	Table1     string `arg:"required,help:table name for db1"`
	Where1     string `arg:"help:SQL filter injected in WHERE clause of db1"`
	Pos1       string `arg:"required,help:reference point for reads of db1; one of 5p or 3p"`
	Collapse1  bool   `arg:"help:collapse reads that have the same pos1"`
	Random1    bool   `arg:"help:randomise pos1"`
	FASTA      string `arg:"help:FASTA file with reference sequences; enables randomization preserving nt content"`
	WinUp      int    `arg:"--win-up,help:max upstream distance from pos to preserve nt content"`
	WinDown    int    `arg:"--win-down,help:max downstream distance from pos to preserve nt content"`
	DownSample bool   `arg:"help:downsample db1 so that content at pos1 is 25% for all nucleotides"`
	DB2        string `arg:"required,help:SQLite3 database 2"`
	Table2     string `arg:"required,help:table name for db2"`
	Where2     string `arg:"help:SQL filter injected in WHERE clause of db2"`
	Pos2       string `arg:"required,help:reference point for reads of db2; one of 5p or 3p"`
	Collapse2  bool   `arg:"help:collapse reads that have the same pos2"`
	BED        string `arg:"required,help:BED file with feats (at least 6 columns)"`
	Span       int    `arg:"required,help:maximum distance of compared pos"`
	Offset     int    `arg:"help:shrink offset 5' and 3' of BED feats"`
	ByFeat     bool   `arg:"--by-feat,help:report distribution by feature"`
	Anti       bool   `arg:"help:Compare reads on opposite instead of same orientation"`
	QArea      int    `arg:"--q-area,help:area (nts) around feats to extend query; use a big number."`
}

// Version returns the program version.
func (Opts) Version() string {
	return "relative-pos-distro 0.7"
}

// Description returns an extended description of the program.
func (Opts) Description() string {
	return "Measure the distribution of relative read positions in database 1 against database 2. Prints the number of read pairs at each relative position, the total number of possible pairs and the total number of reads in each database. Positive relative positions indicate read 1 is downstream of read 2. The counting is done by filtering reads in database 2 within the provided features. Provided SQL filters will apply to all counts."
}

func main() {
	var err error
	var db1, db2 *sqlx.DB

	var opts Opts
	opts.QArea = 1000
	p := arg.MustParse(&opts)
	if opts.Pos1 != "5p" && opts.Pos1 != "3p" {
		p.Fail("--pos1 must be either 5p or 3p")
	}
	if opts.Pos2 != "5p" && opts.Pos2 != "3p" {
		p.Fail("--pos2 must be either 5p or 3p")
	}

	// open database connections.
	if db1, err = sqlx.Connect("sqlite3", opts.DB1); err != nil {
		log.Fatal(err)
	}
	if db2, err = sqlx.Connect("sqlite3", opts.DB2); err != nil {
		log.Fatal(err)
	}

	// create select decorators.
	decors1 := []BuilderDecorator{Table(opts.Table1), Where(opts.Where1)}
	decors2 := []BuilderDecorator{Table(opts.Table2), Where(opts.Where2)}

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
			seqs[s.Name()] = bytes.ToUpper(alphabet.LettersToBytes(s.Seq))
		}
	}

	// goroutine that sends each reference as a job to jobs.
	jobs := make(chan job)
	go func() {
		for bedS.Next() {
			if bedS.Error() != nil {
				log.Fatal("error reading bed: ", bedS.Error())
			}
			b := bedS.Feat()

			seq, ok := seqs[b.Location().Name()]
			if !ok && opts.FASTA != "" {
				continue
			}
			jobs <- job{
				opts:    opts,
				ref:     b,
				seq:     seq,
				db1:     db1,
				db2:     db2,
				decors1: decors1,
				decors2: decors2,
			}
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
	if opts.ByFeat == true {
		fmt.Printf("ref\tpos\tpairs\treadCount1\treadCount2\n")
		for res := range results {
			for i := -opts.Span; i <= opts.Span; i++ {
				fmt.Printf("%s\t%d\t%d\t%d\t%d\n",
					res.job.ref.Location().Name(), i, res.hist[i], res.count1, res.count2)
			}
		}
	} else {
		var totalCount1, totalCount2 int
		aggrHist := make(map[int]uint)
		for res := range results {
			for k, v := range res.hist {
				aggrHist[k] += v
			}
			totalCount1 += res.count1
			totalCount2 += res.count2
		}

		fmt.Printf("pos\tpairs\treadCount1\treadCount2\n")
		for i := -opts.Span; i <= opts.Span; i++ {
			fmt.Printf("%d\t%d\t%d\t%d\n", i, aggrHist[i], totalCount1, totalCount2)
		}
	}
}

func worker(id int, jobs <-chan job, results chan<- result) {
	for job := range jobs {
		var err error
		ref := job.ref

		if ref.Len() < 2*job.opts.Offset {
			continue
		}

		// get position extracting function
		getPos1 := htsdb.Head
		if job.opts.Pos1 == "3p" {
			getPos1 = htsdb.Tail
		}
		getPos2 := htsdb.Head
		if job.opts.Pos2 == "3p" {
			getPos2 = htsdb.Tail
		}

		o, ok := ref.(feat.Orienter)
		if !ok {
			log.Fatal("bed feature does not have defined orientation")
		}

		ori := o.Orientation()
		start := ref.Start() - job.opts.QArea - job.opts.Span
		stop := ref.End() - 1 + job.opts.QArea + job.opts.Span
		rname := ref.Location().Name()
		lowPosLim := ref.Start() - job.opts.Span
		topPosLim := ref.End() + job.opts.Span - 1
		if lowPosLim < 0 {
			lowPosLim = 0
		}
		if job.opts.FASTA != "" && topPosLim >= len(job.seq) {
			topPosLim = len(job.seq) - 1
		}

		// assemble sqlx select builders
		rangeDec := Where(
			"strand = ? AND rname = ? AND START BETWEEN ? AND ? AND STOP BETWEEN ? AND ?")
		readsB1 := DecorateBuilder(htsdb.RangeBuilder, append(job.decors1, rangeDec)...)
		readsB2 := DecorateBuilder(htsdb.RangeBuilder, append(job.decors2, rangeDec)...)

		// prepare statements.
		var readsStmt1, readsStmt2 *sqlx.Stmt
		if readsStmt1, err = prepareStmt(readsB1, job.db1); err != nil {
			log.Fatal(err)
		}
		if readsStmt2, err = prepareStmt(readsB2, job.db2); err != nil {
			log.Fatal(err)
		}

		wig := make([]int, topPosLim+1)
		hist := make(map[int]uint)
		visited := make(map[int]bool)
		var r htsdb.Range
		var count1, count2 int

		// loop on reads in db1.
		ori1 := ori
		if job.opts.Anti == true {
			ori1 = -1 * ori1
		}
		rows1, err := readsStmt1.Queryx(ori1, rname, start, stop, start, stop)
		if err != nil {
			log.Fatal(err)
		}
		for rows1.Next() {
			if err = rows1.StructScan(&r); err != nil {
				log.Fatal(err)
			}
			pos := getPos1(&r, ori1)
			if pos < lowPosLim || pos > topPosLim {
				continue
			}
			if wig[pos] > 0 && job.opts.Collapse1 {
				continue
			}
			count1++
			wig[pos]++
		}

		if job.opts.Random1 {
			if job.opts.FASTA != "" {
				var b = make([]int, len(wig))
				copy(b, wig)
				rand.PermKeepByte(
					wig[lowPosLim:topPosLim],
					job.seq[lowPosLim:topPosLim],
					-job.opts.WinUp, job.opts.WinDown, nil)
			} else {
				rand.Perm(wig[lowPosLim:topPosLim], nil)
			}
		}
		if job.opts.DownSample {
			downSampleKeepByte(wig[lowPosLim:topPosLim], job.seq[lowPosLim:topPosLim])
		}

		// loop on reads in db2.
		rows2, err := readsStmt2.Queryx(ori, rname, start, stop, start, stop)
		if err != nil {
			log.Fatal(err)
		}
		for rows2.Next() {
			if err = rows2.StructScan(&r); err != nil {
				log.Fatal(err)
			}
			pos := getPos2(&r, ori)
			if pos < ref.Start()+job.opts.Offset || pos >= ref.End()-job.opts.Offset {
				continue
			}
			if visited[pos] && job.opts.Collapse2 {
				continue
			}
			visited[pos] = true
			count2++
			for relPos := -job.opts.Span; relPos <= job.opts.Span; relPos++ {
				if pos+relPos < 0 {
					continue
				}
				hist[relPos*int(ori)] += uint(wig[pos+relPos])
			}
		}

		// enqueue in results channel
		results <- result{hist: hist, job: job, count1: count1, count2: count2}
	}
}

func downSampleKeepByte(wig []int, seq []byte) {
	if len(wig) != len(seq) {
		log.Fatal("error: length for wig and seq differ:", len(wig), len(seq))
	}

	var sample = make(map[byte][]int)
	for i, v := range wig {
		if v == 0 {
			continue
		}
		sample[seq[i]] = append(sample[seq[i]], i)
	}

	if len(sample) != 4 {
		for i := range wig {
			wig[i] = 0
		}
		return
	}

	minElements := -1
	var minByte byte
	for k, v := range sample {
		if minElements < 0 || len(v) < minElements {
			minElements = len(v)
			minByte = k
		}
	}

	for k, v := range sample {
		if k == minByte || len(v) <= minElements {
			continue
		}
		rand.PermPos(v, nil)
		for _, index := range v[minElements:len(v)] {
			wig[index] = 0
		}
	}
}

type job struct {
	opts             Opts
	ref              feat.Feature
	seq              []byte
	decors1, decors2 []BuilderDecorator
	db1, db2         *sqlx.DB
}

type result struct {
	hist   map[int]uint
	count1 int
	count2 int
	job    job
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
