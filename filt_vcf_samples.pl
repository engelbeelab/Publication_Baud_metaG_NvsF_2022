#! /usr/bin/perl -w

#Usage: perl filt_vcf_samples.pl firm5_7_corecov_coord.txt freebayes_191105_core.vcf

my %ref_ids = (  
'firm5_1' => 'Ga0133563',
'firm5_2' => 'Ga0133561',
'firm5_3' => 'Ga0133562',
'firm5_4' => 'Ga0133564',
'firm5_bombus_1' => 'Ga0226852',
'firm5_bombus_2' => 'Ga0226847',
'firm5_7' => 'Ga0312303',
'gilli_1' => 'GAPWK',
'gilli_2' => 'Ga0133555',
'gilli_3' => 'Ga0133557',
'gilli_4' => 'Ga0307802',
'gilli_5' => 'Ga0307800',
'gilli_6' => 'Ga0307803',
'gilli_apis_andre' => 'A9G09', #No complete ref
'gilli_apis_dorsa' => 'A9G13', #No complete ref
'gilli_bombus' => 'Ga0227302', #No complete ref
'snod_1' => 'SALWKB2',
'snod_2' => 'BHC53', #No complete ref, contigs have been re-ordered
'snod_bombus' => 'SALWKB12',#No complete ref
'bifido_1.1' => 'H3T92',
'bifido_1.2' => 'H3U89',
'bifido_1.3' => 'BAST',
'bifido_1.4' => 'H3T91',
'bifido_1.5' => 'Ga0133553',
'bifido_1_cerana' => 'Ga0312305',
'bifido_2' => 'BINDI',
'firm4_1' => 'Ga0326456', 
'firm4_2' => 'Ga0072399',
'fper_1' => 'Ga0077910',
'bapis' => 'BBC0122',
'api_1' => 'Ga0312307',
'api_apis_dorsa' => 'C4S76', #No complete ref
'api_bombus' =>  'Ga0061079', #No complete ref
'com_1' => 'D9V35',
'com_drosophila' => 'CIN', #No complete ref
'com_monarch' => 'Ga0248239', #No complete ref
'bom_apis_melli' => 'Ga0216357', #No complete ref
'bom_bombus' => 'Ga0308518', #No complete ref
'bom_1' => 'Ga0372762',
'lkun' => 'Ga0098617',
    );

my $sdp = substr $ARGV[0],0,-18;
my $genome_id = $ref_ids{$sdp};

#Read the *coord.txt file, and get the names of samples with sufficient coverage.
open(FILE, $ARGV[0]) or die
    "Cant open $ARGV[0]!";
my @high_cov_samples;
while(<FILE>) {
    chomp;
    next if ($_ =~ /Cluster/);
    my @split = split("\t",$_);
    my $sample = $split[1];
    my $ter_cov = $split[2];
    if ($ter_cov >= 20) {
	push @high_cov_samples,$sample;
    }
}
close FILE;

#Read the vcf-file, and create a temporary subset vcf-file, containing data for the SDP of interest

open OUTFILE, '>temp.vcf' or die $!;
open(FILE, $ARGV[1]) or die
    "Cant open $ARGV[1]!";
while(<FILE>) {
    chomp;
    if ($_ =~ /#/) {
	print OUTFILE $_,"\n";
    }
    else {
	my @split = split("\t",$_);
	if ($split[0] eq $genome_id) {
	    print OUTFILE $_,"\n";
	}
    }
}  
close FILE;

#Print the list of samples with high coverage, as input for bash-script
print join(" ",@high_cov_samples);

