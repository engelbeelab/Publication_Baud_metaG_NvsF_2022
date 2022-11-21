#! /usr/bin/perl -w
use List::Util qw( min max );
use Pod::Usage;
use Getopt::Long;

pod2usage(-verbose => 2) if ((@ARGV==0)&&(-t STDIN));

=head1 SYNOPSIS

perl filter_sam_aln_length.pl my.sam > my_filt.sam

=head1 DESCRIPTION

This script is used to filter a sam-file according to alignment length only (alignment length > 50). The information is parsed from the CIGAR within the samfile. To use on bam-file, pipe the data with samtools:

samtools view -h my.bam | perl filter_sam_aln_length.pl - > my_filt.sam

=cut

GetOptions(
    'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
) or pod2usage(2);

open (FILE, $ARGV[0]) or
  die "Can't open $ARGV[0]: $!";

my @match_numbers;
while(<FILE>) {
    chomp;
    my @split = split(" ",$_);
    if (/^\@/) {
	print $_,"\n";
    }
    else {
	my $cigar = $split[5];
	my @matches = ($cigar =~ m/([0-9]+M)/g);
	my $nb_matches = @matches;	
	if ($nb_matches > 0) {
	    foreach my $match(@matches) {
		push @match_numbers, substr $match,0,-1;
	    }
	    my $max = max @match_numbers;
	    if ($max <= 50) {
		print $_,"\n";
	    }
	}
	else {
		print $_,"\n";
	}
	@match_numbers = ();
    }
}
