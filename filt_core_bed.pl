#! /usr/bin/perl -w

#Usage: perl filt_core_bed.pl api_1_corecov_coord.txt
#1: read the coord-file, save names of samples with sufficient coverage
#2: read the corecov-file, and check the coverage for all genes for samples with sufficient coverage. Save genes with sufficient coverage.
#3: read the bed-file, and produce a reduced version, using gene-names from step

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
'gilli_apis_andre' => 'A9G03', #No complete ref
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
'bifido_2' => 'BINDI',
'bifido_cerana' => 'Ga0312305',
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

my $sdp = substr $ARGV[0], 0, -18;

#Open the *coord-file, containing the terminus coverage per sample, save names of samples with at least 20x ter-cov
my %high_cov_samples;
open(FILE, $ARGV[0]) or die
    "Cant open $ARGV[0]!";
while(<FILE>) {
    chomp;
    next if ($_ =~ /Cluster/);
    my @split = split("\t",$_);
    my $sample = $split[1];
    my $ter_cov = $split[2];
    if ($ter_cov >= 20) {
	$high_cov_samples{$sample} =$ter_cov;
    }
}
close FILE;

#Open the *corecov-file, containing the gene coverage per core gene, for each sample. For the samples with more than 20x ter-cov, check if the gene coverage is at least 10x. If it is, save the gene-id.
my %low_cov_genes;
my $corecov_file = $sdp."_corecov.txt";
open(FILE, $corecov_file) or die
    "Cant open $corecov_file!";
while(<FILE>) {
    chomp;
    next if ($_ =~ /SDP/);
    my @split = split("\t",$_);
    my $sample = $split[1];
    if (exists $high_cov_samples{$sample}) {
	my $gene_cov = $split[4];
	my $gene_id = $split[2];
	if ($gene_cov < 10) {
	    $low_cov_genes{$gene_id} = 1;
	}
    }
}
close FILE;

#Open the bed-file containing the core-genes, and create a reduced version containing genes with sufficient coverage (i.e. the real core-genes)

my $bed_file = $ref_ids{$sdp}."_core_red.bed";
my $outfile = $sdp."_core_filt.bed";
open OUTFILE, '>',$outfile or die $!;
open(FILE, $bed_file) or die
    "Cant open $bed_file!";
while(<FILE>) {
    chomp;
    my @split = split("\t",$_);
    my $gene_id = $split[3];
    unless (exists $low_cov_genes{$gene_id}) {
	print OUTFILE $_,"\n";
    }
}
close FILE;
close OUTFILE;
