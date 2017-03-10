#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 2;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $gtf, @bams, @options, $outfile, );

my $exec = $plugin_path . "/../bin/quantifyReadLengthOverGtf.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/quantifyReadLengthOverGtf";
unless ( -d "$plugin_path/data/output/") {
	mkdir ( "$plugin_path/data/output/" );
}
if ( -d $outpath ) {
	system("rm -Rf $outpath/*");
}

$outfile= "$outpath/run1";

@bams = map { "$plugin_path/data/$_.bam" } map { $_."1", $_."2", $_."3"} 'A', 'B';

foreach ( @bams ) {
	ok ( -f $_, "bam file $_ exists" );
}

$gtf =  "$plugin_path/data/Small.gtf";

ok ( -f $gtf, "gtf file $gtf exists" );

@options = ( 'quantify_on', 'gene' );

my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -gtf " . $gtf 
. " -bams " . join(' ', @bams )
. " -options " . join(' ', @options )
. " -outfile " . $outfile 
#. " -debug"
;
my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";

#print "\$exp = ".root->print_perl_var_def($value ).";\n";