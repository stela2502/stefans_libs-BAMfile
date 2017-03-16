#! /usr/bin/perl -w

=head1 LICENCE

  Copyright (C) 2017-03-06 Stefan Lang

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

    quantifyReadLengthOverGtf.pl
       -bams     :a list of bam files you want to get analyzed
       -options  :  
 dist_from_feature_start:how far from the feature start may the read start (default 1bp)
        size_fractions: e.g. "20 40" (default) or "15 30 50" 
                        to report the reads in lengt of below 20 below 40 or above 40
                        or below 15, below 30, below 50 or more than 50
           quantify_on: which feature to select for quantification from the gtf file (default 'exon')
             report_on: which column in the gtf file object to report on (default 'gene_id')
          slice_length: change the sice of the feature matching objects (default 5e+6)      

       -gtf       :the gtf file with the features to quantify
       -outfile   :the outfile (a tab separated table)


       -help           :print this help
       -debug          :verbose output
   
=head1 DESCRIPTION

  This tool creates a histogram over the read length on a spcific type of entries in the gtf file. It has been implemented to identify whether tRNA has been fragmented or not and whether these fragments start at position 0 in the RNA or not.

  To get further help use 'quantifyReadLengthOverGtf.pl -help' at the comman line.

=cut

use Getopt::Long;
use Pod::Usage;

use stefans_libs::BAMfile;
use stefans_libs::file_readers::gtf_file;

use strict;
use warnings;

use FindBin;
my $plugin_path = "$FindBin::Bin";

my $VERSION = 'v1.0';

my ( $help, $debug, $database, @bams, $options, @options, $gtf, $outfile );

Getopt::Long::GetOptions(
	"-bams=s{,}"    => \@bams,
	"-options=s{,}" => \@options,
	"-gtf=s"        => \$gtf,
	"-outfile=s"    => \$outfile,

	"-help"  => \$help,
	"-debug" => \$debug
);

my $warn  = '';
my $error = '';

unless ( defined $bams[0] ) {
	$error .= "the cmd line switch -bams is undefined!\n";
}
unless ( defined $options[0] ) {
	$warn .= "the cmd line switch -options is undefined!\n";
}
unless ( defined $gtf ) {
	$error .= "the cmd line switch -gtf is undefined!\n";
}
unless ( defined $outfile ) {
	$error .= "the cmd line switch -outfile is undefined!\n";
}

if ($help) {
	print helpString();
	exit;
}

if ( $error =~ m/\w/ ) {
	helpString($error);
	exit;
}

sub helpString {
	my $errorMessage = shift;
	$errorMessage = ' ' unless ( defined $errorMessage );
	print "$errorMessage.\n";
	pod2usage( q(-verbose) => 1 );
}

### initialize default options:

#$options->{'n'} ||= 10;

###

my ($task_description);

$task_description .= 'perl ' . $plugin_path . '/quantifyReadLengthOverGtf.pl';
$task_description .= ' -bams "' . join( '" "', @bams ) . '"'
  if ( defined $bams[0] );
$task_description .= ' -options "' . join( '" "', @options ) . '"'
  if ( defined $options[0] );
$task_description .= " -gtf '$gtf'"         if ( defined $gtf );
$task_description .= " -outfile '$outfile'" if ( defined $outfile );

for ( my $i = 0 ; $i < @options ; $i += 2 ) {
	$options[ $i + 1 ] =~ s/\n/ /g;
	$options->{ $options[$i] } = $options[ $i + 1 ];
}
###### default options ########
$options->{'quantify_on'}    ||= 'exon';
$options->{'report_on'}      ||= 'gene_id';
$options->{'size_fractions'} ||= "20 40";
$options->{'size_fractions'} = [ split( " ", $options->{'size_fractions'} ) ];
$options->{'dist_from_feature_start'} = 1;
$options->{'slice_length'} ||= 5e+6;
##############################
my $fm = root->filemap($outfile);
mkdir( $fm->{'path'} ) unless ( -d $fm->{'path'} );

open( LOG, ">$outfile.log" )
  or die "I could not open the log file '$outfile.log'\n" . $!;
print LOG $task_description . "\n";
close(LOG);

my (
	@matching_Hashes, $runs,         $sample_table, $sample_row,
	$result,          @matching_IDs, $N,            $Seq,
	$sample_name,     $self
);
$runs         = 0;
$self         = {'chr' => ''};                  # the local script variables stash
$sample_table = data_table->new();
$sample_table->add_column( 'filename' => @bams );
$sample_table->Add_2_Header( [ 'mapped', 'in gene', 'unmapped' ] );

$result = data_table->new();
$result->Add_2_Header( ['Gene_ID'] );

my $gtf_obj = stefans_libs::file_readers::gtf_file->new();

$gtf_obj->read_file($gtf);

my $quantifer =
  $gtf_obj->select_where( 'feature',
	sub { $_[0] eq $options->{'quantify_on'} } );
$quantifer->{'slice_length'} = $options->{'slice_length'};

my @bam_line;

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

	unless ( $bam_line[2] =~ m/^chr/ ) {
		@{ @{ $sample_table->{'data'} }[ $BAM_file->{'tmp_sample_row'} ] }[3]++;
		return;
	}
	else {
		@{ @{ $sample_table->{'data'} }[ $BAM_file->{'tmp_sample_row'} ] }[1]++;
	}

	## start with the real matching
	@matching_IDs = &get_matching_ids( $quantifer, $bam_line[2], $bam_line[3] );

	if ( @matching_IDs == 0 ) {    ## no match to any gene / exon
		    #print "BAM line did not match to features - next\n";
		return;
	}
	@{ @{ $sample_table->{'data'} }[$sample_row] }[2]++;

	@matching_Hashes = &get_reporter_hashes( $quantifer, @matching_IDs );

	## now I need to check if the read starts at the start of the feature
	my @tmp;
	foreach my $gtf_hash (@matching_Hashes) {
		$N   = 0;
		$Seq = 0;
		map { $N   += $_ } $bam_line[5] =~ m/(\d+)N/g;
		map { $Seq += $_ } $bam_line[5] =~ m/(\d+)M/g;

#		warn
#"$Seq,  $bam_line[2]:$gtf_hash->{'strand'}: abs($bam_line[3] - $gtf_hash->{start}) "
#		  . abs( $bam_line[3] - $gtf_hash->{start} )
#		  . " <= $options->{'dist_from_feature_start'}\n";
#		warn
#"$Seq,  $bam_line[2]:$gtf_hash->{'strand'}: abs(($bam_line[3]+ $Seq) - $gtf_hash->{end}) "
#		  . abs( ( $bam_line[3] + $Seq ) - $gtf_hash->{end} )
#		  . " <= $options->{'dist_from_feature_start'}\n";
		if ( $N > $Seq ) {
				warn "this can not be used here :-(:\n"
				  . join( "\t", @bam_line ) . "\n";
				next;    ## this can not be used here :-(
		}
		if ( $gtf_hash->{'strand'} eq "+" ) {
			if (
				abs( $bam_line[3] - $gtf_hash->{start} ) <=
				$options->{'dist_from_feature_start'} )
			{
				&add_to_summary( $BAM_file->{'tmp_sample_name'},
					$gtf_hash->{ $options->{'report_on'} }, $Seq );
			}
			elsif (
				abs( ( $bam_line[3] + $Seq ) - $gtf_hash->{end} ) <=
				$options->{'dist_from_feature_start'} )
			{    ## close to end
				&add_to_summary( "$BAM_file->{'tmp_sample_name'} from end",
					$gtf_hash->{ $options->{'report_on'} }, $Seq );
			}
		}
		else {
			if (
				abs( ( $bam_line[3] + $Seq ) - $gtf_hash->{end} ) <=
				$options->{'dist_from_feature_start'} )
			{
				&add_to_summary( $BAM_file->{'tmp_sample_name'},
					$gtf_hash->{ $options->{'report_on'} }, $Seq );
			}
			elsif (
				abs( $bam_line[3] - $gtf_hash->{start} ) <=
				$options->{'dist_from_feature_start'} )
			{            ## colse to end
				&add_to_summary( "$BAM_file->{'tmp_sample_name'} from end",
					$gtf_hash->{ $options->{'report_on'} }, $Seq );
			}
		}
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

sub get_matching_ids {
	my ( $gtf, $chr, $start ) = @_;
	my ( @IDS, $nextID );
	unless ( $self->{'chr'} eq $chr ) {
		$self->{'chr'} = $chr;
		$self->{'end'} = 0;
	}

	if ( $self->{'end'} <= $start ) {
		@IDS = $gtf->efficient_match_chr_position_plus_one( $chr, $start );

		return
		  if ( @IDS == 0 or !defined( $IDS[0] ) )
		  ;    ## pdl has no more matching entries
		$nextID = pop(@IDS);
		unless ( defined $nextID )
		{      ## the end of the chromosome annotation data has been reached
			$self->{'next_start'} =
			  ( $gtf->get_chr_subID_4_start($start) + 1 ) *
			  $gtf->{'slice_length'};
			return;
		}
		$self->{'next_start'} =
		  @{ @{ $gtf->{'data'} }[$nextID] }[ $gtf->Header_Position('start') ];
		$self->{'end'} = &lmin( &get_chr_end_4_ids( $gtf, @IDS ) );
		$self->{'last_IDS'} = \@IDS;
	}
	elsif ( $start < $self->{'next_start'} ) {
		@IDS = ();    ## no match needed...
	}
	else {
		@IDS = @{ $self->{'last_IDS'} };
	}
	return @IDS;
}

sub get_reporter_hashes {
	my ( $gtf, @lines ) = @_;
	return &unique_hash( map { $gtf->get_line_asHash($_) } @lines );
}

sub unique_hash {
	my $h;
	map {
		Carp::confess(
			    "missing $options->{'report_on'} in one \@_ entry\$exp = "
			  . root->print_perl_var_def($_)
			  . ";\n" )
		  unless ( defined $_->{ $options->{'report_on'} } );
		$h->{ $_->{ $options->{'report_on'} } } = $_
	} @_;
	return values %$h;
}

sub get_reporter_ids {
	my ( $gtf, @lines ) = @_;
	return &unique(
		map {
			@{ @{ $gtf->{'data'} }[$_] }
			  [ $gtf->Header_Position( $options->{'report_on'} ) ];
		} @lines
	);
}

sub unique {
	my $h;
	map { $h->{$_}++ } @_;
	return sort keys %$h;
}

sub lmax {
	return 0 if ( @_ == 0 );
	my $max = shift;
	foreach (@_) {
		$_ ||= 0;
		$max = $_ if ( $_ > $max );
	}
	return $max;
}

sub lmin {
	return 0 if ( @_ == 0 );
	my $min = shift;

	#Carp::confess ( 'min: '.join(", ",@_)."\n");
	foreach (@_) {
		$min = $_ if ( $_ < $min );
	}
	return $min;
}

sub get_chr_end_4_ids {
	my ( $gtf, @lines ) = @_;
	return
	  map { @{ @{ $gtf->{'data'} }[$_] }[ $gtf->Header_Position('end') ]; }
	  @lines;
}
