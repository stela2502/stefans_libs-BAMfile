package stefans_libs::BAMfile;

#use FindBin;
#use lib "$FindBin::Bin/../lib/";
use strict;
use warnings;

use stefans_libs::root;
use stefans_libs::file_readers::bed_file;
=head1 LICENCE

  Copyright (C) 2016-10-13 Stefan Lang

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


=for comment

This document is in Pod format.  To read this, use a Pod formatter,
like 'perldoc perlpod'.

=head1 NAME

stefans_libs::BAMfile

=head1 DESCRIPTION

SIMPLE bam file interface that uses samtools view in the back - not a compressed file reader.

=head2 depends on


=cut


=head1 METHODS

=head2 new ( $hash )

new returns a new object reference of the class stefans_libs::BAMfile.
All entries of the hash will be copied into the objects hash - be careful t use that right!

=cut

sub new{

	my ( $class, $hash ) = @_;

	my ( $self );

	$self = {
  	};
  	foreach ( keys %{$hash} ) {
  		$self-> {$_} = $hash->{$_};
  	}

  	bless $self, $class  if ( $class eq "stefans_libs::BAMfile" );

  	return $self;

}


sub open_file {
	my ( $self, $file ) = @_;
	Carp::confess ( "Not a file '$file'\n$!\n" ) unless ( -f $file ) ;
	my $F;
	if ( $file =~ m/bam$/ ) {
		open ( $F, "samtools view -h $file |") ; 
	}elsif ( $file =~ m/sam$/ ){
		open ( $F, "<$file" );
	}
	else {
		Carp::confess ( "Sorry, but ".ref($self)." can not open this file for you : $file\nFormat not supported (not *.bam or *.sam)\n");
	}
	return $F;
}

=head2 dropUMI_from_bam ( sortedBamFile, outfile )

=cut
sub dropUMI_from_bam {
	my ( $self, $sortedBamFile, $outfile ) =@_;
	my $F = $self->open_file( $sortedBamFile );
	my $OUT;
	$outfile .= ".bam" unless ( $outfile =~ m/.bam$/);
	open ( $OUT, "| samtools view -Sb - > $outfile") or die $!;
	my ($entries, $old_entries, $last_position, $last_chr, @line,@tmp, $UMI);
	$last_position = -100;
	$last_chr = 'none';
	$entries = {};
	$self->{dropped} = {};
	while( <$F> ) {
		if ( $_ =~ m/^\@/ ){ ## process header
			print $OUT $_;
		}else {
			#M04223:22:000000000-AUY26:1:1101:6136:18821:ACAAC	272	chr1	629499	0	4S38M	*	0	0	TTGGCTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCA	HHHHHHHHHHHHHHHHHHHHGHGHHGGHGGGGGGGGGGGFFF	AS:i:-4	ZS:i:-4	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:38	YT:Z:UU	NH:i:2
			
			@line=split( "\t", $_);
			next if ( $line[2] eq "*" ); ## kick not mapped
			## check that the position has not changed (+-5bp?)
			if ( ($last_position < ($line[3] - 5)) or !($last_chr eq $line[2]) ) {
				$last_position = $line[3];
				$last_chr = $line[2];
				$old_entries = $entries; 
				$entries = {};
				$old_entries = {} if ( $last_position < ($line[3] - 10));
			}
			@tmp = split( ":", $line[0]);
			$UMI = pop (@tmp);
			if ($old_entries->{$UMI} or $entries->{$UMI} ){
				$self->{'dropped'} ->{$UMI} ||=0;
				$self->{'dropped'} ->{$UMI} ++;
				next;
			}
			$entries->{$UMI} = 1;
			print $OUT $_;
		}
		
	}
	close ( $F);
	close ( $OUT );
	return $self;
	
}

sub filter_file {
	my ( $self, $fname, $filter ) = @_;
	my $file = $self->open_file($fname);
	while ( <$file> ) {
		chomp($_);
		&{$filter}( $self, $_ );
	}
	close ( $file);
	return $self;
}

sub select_4_str {
	my ( $self, $file, $str, $where, $outfile) =@_;
	$where ||=  0;
	my $OUT;
	open ( $OUT, "| samtools view -Sb - > $outfile") or die "I could not create the outfile $outfile\n $!\n";
	#open ( $OUT,">$outfile" );
	my $function = sub {
		my @tmp;
		if ( $_ =~ m/^\@/ ){ ## process header
			print $OUT $_."\n";
		}else {
			@tmp=split( "\t", $_);
			if ( $tmp[$where] =~ m/$str/ ){
				$self->{'OK'} ++;
				print $OUT join("\t", @tmp)."\n";
			}else {
				$self->{'filtered'} ++;
			}
		}
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file($file, $function );
	close ( $OUT );
	print "$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}

=head2 search_flag_binary( $file, [0,0,1,0,1], $outfile)

Searches if the 1 bits are set in the flag.
If any bit is not set, the hit will be rejected.
0 bits will not be checked.

=cut

sub search_flag_binary {
	my ( $self, $file, $bits, $outfile) =@_;
	$bits ||=  [0,0,0,1]; ## reverse complement 
	my $OUT;
	open ( $OUT, "| samtools view -Sb - > $outfile")  or die "I could not create the outfile $outfile\n $!\n";
	my $function = sub {
		my @tmp;
		if ( $_ =~ m/^\@/ ){ ## process header
			print $OUT $_."\n";
		}else {
			@tmp=split( "\t", $_);
			my @b = split(//, sprintf ( "%b", $tmp[1] ));
			my $ok = 1;
			foreach (my $i = 0; $i < @$bits; $i++ ) {
				if ( @$bits[$i] ) {
					$ok = 0 unless ( $b[$i] );
				}
			}
			if ( ! $ok ) { 
				$self->{'filtered'} ++;
			}
			else {
				$self->{'OK'} ++;
				print $OUT join("\t", @tmp)."\n";
			}
		}
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file($file, $function );
	close ( $OUT );
	print "$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}

sub select_single_match {
	my ( $self, $file, $max_hits, $outfile) =@_;
	my $OUT;
	$max_hits ||= 1;
	open ( $OUT, "| samtools view -Sb - > $outfile")  or die "I could not create the outfile $outfile\n $!\n";
	my $function = sub {
		my @tmp;
		if ( $_ =~ m/^\@/ ){ ## process header
			print $OUT $_."\n";
		}else {
			if ( $_ =~ m/NH:i:(\d+)/ and $1 > $max_hits ) { 
				$self->{'filtered'} ++;
			}
			else {
				$self->{'OK'} ++;
				print $OUT $_."\n";
			}
		}
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file($file, $function );
	close ( $OUT );
	print "$self->{'OK'} OK fastq entries, $self->{'filtered'} filtered.\n";
	return $self;
}

sub toBed {
	my ( $self, $file, $outfile ) = @_;
	my $ret = stefans_libs::file_readers::bed_file->new();
	#Plate_1_G04_primer1	0	chr1	4493818	60	20M	*	0	0	GGGGAAATAGGAAGGCTGAA	IIIIIIIIIIIIIIIIIIII	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:20	YT:Z:UU	NH:i:1
	
	my $function = sub {
		my @line;
		my $add;
		if ( $_ =~ m/^\@/ ){ ## process header
		}else {
			my @line= split( "\t",$_);
			$add = 0;
			map { $add += $_ } split(/[MIDNSHPX\=*]/, $line[5] ); 
			push( @{$ret->{'data'}}, [$line[2],$line[3],$line[3]+ $add, $line[0]] );
		}
	};
	$self->{'OK'} = $self->{'filtered'} = 0;
	$self->filter_file($file, $function );
	if ( defined $outfile ) {
		$ret->write_file($outfile);
	}
	return $ret;
}


1;
