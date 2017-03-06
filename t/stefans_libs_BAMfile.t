#! /usr/bin/perl
system( 'module load SAMtools/1.3.1-HTSlib-1.3.1');
use strict;
use warnings;
use Test::More tests => 9;
BEGIN { use_ok 'stefans_libs::BAMfile' }
use stefans_libs::root;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my ( $value, @values, $exp, @tmp );
my $OBJ = stefans_libs::BAMfile -> new({'debug' => 1});
is_deeply ( ref($OBJ) , 'stefans_libs::BAMfile', 'simple test of function stefans_libs::BAMfile -> new() ');

$value = $OBJ -> open_file ($plugin_path."/data/problems.bam" );
is_deeply ( ref($value), 'GLOB', 'opened bam file' );
@values = <$value>;
close ( $value );
$value = [shift(@values), pop(@values )];
#print "\$exp = ".root->print_perl_var_def($value ).";\n";
$exp = [ '@HD	VN:1.0	SO:coordinate
', 'M04223:22:000000000-AUY26:1:2114:13545:28679:ACACC	4	*	0	0	*	*	0	0	TCCTCGTTTGTATAGTGGTGAGTATCCCCGACGGGGAGCCAA	FFBGGGGGGGGGGHHHHGHFFHHHHHHHGGEGGGGGGGGGEG	YT:Z:UU
' ];

is_deeply( $value, $exp, "bam file content" );


$value = $OBJ -> open_file ($plugin_path."/data/problems.sam" );
is_deeply ( ref($value), 'GLOB', 'opened sam file' );
@values = <$value>;
is_deeply( scalar(@values), 9658, "right number of reads + header" );
close ( $value );
$value = [shift(@values), pop(@values )];
is_deeply( $value, $exp, "sam file content" );


$OBJ->dropUMI_from_bam (  $plugin_path."/data/problems.sam", $plugin_path."/data/outpath/dropped_UMIs" );
$value = $OBJ -> open_file ($plugin_path."/data/outpath/dropped_UMIs.bam" );
@values = <$value>;
is_deeply( scalar(@values), 1143, "right number of reads -UMI + header" );
close ( $value );
foreach (@values) {
	next if ( $_ =~ m/^@/ );
	chomp($_);
	push(@tmp, $_);
}
#print "\$exp = ".root->print_perl_var_def( [scalar(@tmp)] ).";\n";

is_deeply( scalar(@tmp), 686, "dropped the right number of duplicates" );


#print "\$exp = ".root->print_perl_var_def($value ).";\n";


