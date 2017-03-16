#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-03-15 Stefan Lang

  This program is free software; you can redistribute it 
  and/or modify it under the terms of the GNU General Public License 
  as published by the Free Software Foundation; 
  either version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License 
  along with this program; if not, see <http://www.gnu.org/licenses/>.

=head1  SYNOPSIS

    quantifyReadLengthOverChr.pl
       -bams     :<please add some info!> you can specify more entries to that
       -options     :<please add some info!> you can specify more entries to that
                         format: key_1 value_1 key_2 value_2 ... key_n value_n
       -outfile       :<please add some info!>


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  more simple version of quantifyReadLengthOverGtf.pl where only the reads mapping to a specific position on each chromosome are counted.

  To get further help use 'quantifyReadLengthOverChr.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::BAMfile;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';


my ( $help, $debug, $database, @bams, $options, @options, $outfile);

Getopt::Long::GetOptions(
       "-bams=s{,}"    => \@bams,
       "-options=s{,}"    => \@options,
	 "-outfile=s"    => \$outfile,

	 "-help"             => \$help,
	 "-debug"            => \$debug
);

my $warn = '';
my $error = '';

unless ( defined $bams[0]) {
	$error .= "the cmd line switch -bams is undefined!\n";
}
unless ( defined $options[0]) {
	$warn .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $outfile) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}


if ( $help ){
	print helpString( ) ;
	exit;
}

if ( $error =~ m/\w/ ){
	helpString($error ) ;
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage); 
	print "$errorMessage.\n";
	pod2usage(q(-verbose) => 1);
}

### initialize default options:

#$options->{'n'} ||= 10;

###


my ( $task_description);

$task_description .= 'perl '.$plugin_path .'/quantifyReadLengthOverChr.pl';
$task_description .= ' -bams "'.join( '" "', @bams ).'"' if ( defined $bams[0]);
$task_description .= ' -options "'.join( '" "', @options ).'"' if ( defined $options[0]);
$task_description .= " -outfile '$outfile'" if (defined $outfile);


for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
$options->{'chr_position'} ||= 1;
$options->{'chr_max_dist'} ||= 1;
$options->{'size_fractions'} ||= "20 40";
$options->{'size_fractions'} = [ split( " ", $options->{'size_fractions'} ) ];
#$options->{'something'} ||= 'default value';
##############################
my $fm = root->filemap( $outfile );
mkdir( $fm->{'path'}) unless ( -d $fm->{'path'} );

open ( LOG , ">$outfile.log") or die $!;
print LOG $task_description."\n";
close ( LOG );


## Do whatever you want!


my (@bam_line,$sample_name, $runs, $sample_table, $N, $Seq, $result, $sample_row) ;
$runs=0;

$result = data_table->new();
$result->Add_2_Header('Gene_ID');
$sample_table = data_table->new();
$sample_table->add_column( 'filename' => @bams );
$sample_table->Add_2_Header( [ 'mapped', 'in gene', 'unmapped' ] );

sub filter {
	my $BAM_file    = shift;
	my $sample_name = $BAM_file->{'tmp_sample_name'};
	my $sample_row  = $BAM_file->{'tmp_sample_row'};

  #Carp::confess ( "I hope the sample name '$sample_name' is a bam filename\n");
	@bam_line = split( "\t", shift );
	return if ( $bam_line[0] =~ m/^\@/ );
	$runs++;

	if ( $runs % 1e4 == 0 ) {
		print ".";
	}

	unless ( $bam_line[2] =~ m/chr/ ) {
		@{ @{ $sample_table->{'data'} }[ $BAM_file->{'tmp_sample_row'} ] }[3]++;
		return;
	}
	else {
		@{ @{ $sample_table->{'data'} }[ $BAM_file->{'tmp_sample_row'} ] }[1]++;
	}


	## now I need to check if the read starts at the start of the feature
	my @tmp;
	
	$N   = 0;
	$Seq = 0;
	map { $N   += $_ } $bam_line[5] =~ m/(\d+)N/g;
	map { $Seq += $_ } $bam_line[5] =~ m/(\d+)M/g;


	if ( $N > $Seq ) {
			warn "this can not be used here :-(:\n"
			  . join( "\t", @bam_line ) . "\n";
			next;    ## this can not be used here :-(
	}
	if ( abs( $bam_line[3] - $options->{'chr_position'} ) <= $options->{'chr_max_dist'} ){
		&add_to_summary( $BAM_file->{'tmp_sample_name'}, $bam_line[2], $Seq );
	}
		

}

my $bam_file = stefans_libs::BAMfile->new();

$| = 1;## turn on autoflush for the process bar

foreach $sample_name (@bams) {
	$result->Add_2_Header(
		[
			(map { "$sample_name $_" } @{ $options->{'size_fractions'} },
			'larger'),
			(map { "$sample_name from end $_" }
			  @{ $options->{'size_fractions'} },
			'larger')
		]
	);
	$sample_row = undef;
	($sample_row) =
	  $sample_table->get_rowNumbers_4_columnName_and_Entry( 'filename',
		$sample_name );
	$bam_file->{'tmp_sample_name'} = $sample_name;
	$bam_file->{'tmp_sample_row'}  = $sample_row;
	$bam_file->filter_file( $sample_name, \&filter );
	print "\ndone with file $sample_name\n";
}
$| = 0;## turn off autoflush

$result->write_file($outfile);
$sample_table->write_file(
	$fm->{'path'} ."/". $fm->{'filename_base'} . '_sampleInfo' );


sub add_to_summary {
	my ( $sampleID, $geneIDs, $read_length ) = @_;

	#warn "add_to_summary: $sampleID $read_length\n";
	my @colnames = map { "$sampleID $_" } @{ $options->{'size_fractions'} },
	  'larger';
	my @col_number = $result->Header_Position( \@colnames );
	map {
		Carp::confess("the col_number for sample $colnames[$_] is undefined\n")
		  unless ( defined $col_number[$_] )
	} 0 .. ( @col_number - 1 );
	$geneIDs = [$geneIDs] unless ( ref($geneIDs) eq "ARRAY" );
	my ( @row_numbers, $row, $added );
	foreach my $gene_id (@$geneIDs) {
		@row_numbers =
		  $result->get_rowNumbers_4_columnName_and_Entry( 'Gene_ID', $gene_id );
		if ( @row_numbers == 0 ) {    ## gene never occured
			    #print "new gene $gene_id detected for sample $sampleID\n";
			my $hash = { 'Gene_ID' => $gene_id };
			$result->AddDataset($hash);
		}
		foreach my $row (@row_numbers) {
			$added = 0;
			for ( my $col_id = 0 ; $col_id < @col_number - 1 ; $col_id++ ) {
				if ( $read_length < @{ $options->{'size_fractions'} }[$col_id] )
				{
#				warn join( "\t", @bam_line )
#				  . "\nI can add size @{$options->{'size_fractions'}}[$col_id] for gene $gene_id and read_length $read_length\n";
					@{ @{ $result->{'data'} }[$row] }[ $col_number[$col_id] ]++;
					$added = 1;
					last;
				}
			}
			unless ($added) {

#			warn join( "\t", @bam_line )
#			  . "\nI can add size 'larger' for gene $gene_id and read_length $read_length\n";
				@{ @{ $result->{'data'} }[$row] }[ $col_number[-1] ]++;
			}
		}
	}

}

