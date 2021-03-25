#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;
use GenOO::TranscriptCollection::Factory;
use PDL::Lite; $PDL::BIGPDL = 0; $PDL::BIGPDL++; # enable huge pdls

# Define and read command line options
my ($opt, $usage) = describe_options(
	'%c %o',
	['gtf=s', 'GTF file with transcripts', { required => 1}],
	['cons-dir=s', 'directory with phastCons files', { required => 1}],
	['cons-file-suffix=s', 'files should be in the form rname+sufix',
		{ default => '.phastCons100way.wigFix.gz'}],
	['chr-pref', 'add "chr" prefix to coords to match phastCons file names'],
	['help|h', 'print usage message and exit'],
);
print($usage->text), exit if $opt->help;

warn "Reading transcripts\n";
my $tr_collection = GenOO::TranscriptCollection::Factory->create('GTF', {
	file => $opt->gtf
})->read_collection;
my @transcripts = sort {$a->rname cmp $b->rname} $tr_collection->all_records;

my %out;
my $open_rname_file = "";
my $rname_pdl;
my $chrpref = '';
if ($opt->chr_pref) {
	$chrpref = "chr";
}
foreach my $tr (@transcripts) {
	my @tr_vals;
	my $rname_file = $opt->cons_dir . $chrpref . $tr->rname . $opt->cons_file_suffix;
	if ($rname_file ne $open_rname_file) {
		warn "Opening file $rname_file\n";
		if (!-e $rname_file) {
			warn "skipping $rname_file\n" if ! -e $rname_file;
			next;
		}
		$rname_pdl = phylop_pdl_for($rname_file);
		$open_rname_file = $rname_file;
	}

	my @exons = @{$tr->exons};
	foreach my $e (@exons) {
		for (my $pos = $e->start; $pos <= $e->stop; $pos++) {
			push @tr_vals, $rname_pdl->at($pos); 
		}
	}
	
	if ($tr->strand == -1) {
		@tr_vals = reverse(@tr_vals);
	}

	my $expectedLen = $tr->exonic_length;
	die if @tr_vals != $expectedLen;

	$out{$tr->id} = '"'.$tr->id.'"' . " : [".join(",", @tr_vals)."]";
}

print "{\n";
print join (",\n", (values %out));
print "\n}\n";

sub phylop_pdl_for {
	my ($rname_file) = @_;

	my $pdl = PDL->zeros(PDL::short(), 300000000);

	open (my $H, "gzip -dc $rname_file |");

	my ($start, $step);
	while (my $line = $H->getline) {
		chomp $line;
		#fixedStep chrom=chr1 start=10918 step=1
		if ($line =~ /^fixedStep\schrom=(.+)\sstart=(.+)\sstep=(.+)/) {
			($start, $step) = ($2 - 1, $3);
		}
		else {
			if ($step == 1) {
				$pdl->set($start, int(1000 * $line));
			} else {
				my $pdl_region = $pdl->slice([$start, $start + $step - 1]) .= int(1000 * $line);
			}
			$start += $step;
		}
	}
	close $H;

	return $pdl;
}
