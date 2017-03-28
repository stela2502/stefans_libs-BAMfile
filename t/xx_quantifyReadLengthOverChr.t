#! /usr/bin/perl
use strict;
use warnings;
use stefans_libs::root;
use Test::More tests => 6;
use stefans_libs::flexible_data_structures::data_table;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, $outfile, @bams, @options, );

my $exec = $plugin_path . "/../bin/quantifyReadLengthOverChr.pl";
ok( -f $exec, 'the script has been found' );
my $outpath = "$plugin_path/data/output/quantifyReadLengthOverChr";
if ( -d $outpath ) {
	system("rm -Rf $outpath");
}

foreach ( 'Simple_Chr_scan2.bam', 'Simple_Chr_scan.bam'){
	push( @bams, $plugin_path."/data/$_" );
	ok (-f $plugin_path."/data/$_", "infile $plugin_path/data/$_ exists");
}

$outfile = "$outpath/run1";


my $cmd =
    "perl -I $plugin_path/../lib  $exec "
. " -outfile " . $outfile 
. " -bams " . join(' ', @bams )
. " -debug";

$cmd .= " -options " . join(' ', @options ) if ( $options[0] );

my $start = time;
system( $cmd );
my $duration = time - $start;
print "Execution time: $duration s\n";

foreach ( 'run1_sampleInfo.xls', 'run1.log', 'run1.xls' ) {
	ok ( -f "$outpath/$_" , "outfile $outpath/$_");
}

@bams = ( $plugin_path."/data/Homo_sapiens_nm-tRNA-Tyr-GTA-chr21-2.bam");
ok ( -f $bams[0], "defined bam file" );



#print "\$exp = ".root->print_perl_var_def($value ).";\n";