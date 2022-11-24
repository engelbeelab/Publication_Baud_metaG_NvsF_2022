#! /usr/bin/perl -w
use Getopt::Long;

#Usage: perl core_cov_red.pl --bam-list list_of_bamfiles.txt

#This script will calculate the core gene coverage when mapping against the reduced database. Bamfiles are indicated with the --bam-list argument, one line per file-name, and these must be sorted, indexed and present in the run directory. Moreover, bed-files containing the filtered core genes for the genomes in the reduced db must be present. Files with core-coverage in longformat are automatically generated ([SDP]_corecov.txt).

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
'gilli_apis_andre' => 'A9G03',
'gilli_apis_dorsa' => 'A9G13',
'gilli_bombus' => 'Ga0227302',
'snod_1' => 'SALWKB2',
'snod_2' => 'BHC53',
'snod_bombus' => 'SALWKB12',
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
'api_apis_dorsa' => 'C4S76',
'api_bombus' => 'Ga0061079',
'com_1' => 'D9V35',
'com_drosophila' => 'CIN',
'com_monarch' => 'Ga0248239',  
'lkun' => 'Ga0098617',
'bom_1' => 'Ga0372762',
'bom_apis_melli' => 'Ga0216357', 
'bom_bombus' => 'Ga0308518', 
    );
my %ref_ids_rv = reverse %ref_ids;

my $bam_list;
GetOptions (
         'b|bam-list:s' => \$bam_list,
         'help' => sub { pod2usage( -exitstatus => 0, -verbose => 2 ) },
    ) or pod2usage(2);
open(FILE, $bam_list) or die
    "Cant open $bam_list!";
my @bam_files;
while(<FILE>) {
    chomp;
    push @bam_files,$_;
}
close FILE;


#Calculate coverage on all core genes in all samples
my %gene_cov;
my %gene_pos;
my %samples;
foreach my $bam_file(@bam_files) {
    my @split_name = split("_",$bam_file);
    my $sample = $split_name[0];
    $samples{$sample}=1;
    print "Processing bam-file: ",$bam_file,"\n";
    foreach my $genome( sort keys %ref_ids_rv) {
	print "\t",$ref_ids_rv{$genome},"\n";
	my $bedfile = $genome."_core_red.bed";
	my @cov_lines = `samtools bedcov $bedfile $bam_file`;
	foreach my $line(@cov_lines) {
	    chomp($line);
	    my @split_line = split("\t",$line);
	    my $gene_start = $split_line[1];
	    my $gene_id = $split_line[3];
	    my $gene_length = $split_line[2]-$split_line[1];
	    my $gene_cov = sprintf("%.2f",($split_line[4]/$gene_length));
	    $gene_cov{$genome}{$gene_id}{$sample} = $gene_cov;
	    $gene_pos{$genome}{$gene_id} = $gene_start;
	}
    }
}
my @samples = sort keys %samples;

#Generate output files in long-format, one file per ref-genome

foreach my $genome( sort keys %ref_ids_rv) {
    my $outfile = $ref_ids_rv{$genome}."_corecov.txt";
    open OUTFILE, '>'.$outfile or die $!;
    print OUTFILE "SDP\tSample\tGene_id\tRef_pos\tCoverage\n";
    my $sdp = $ref_ids_rv{$genome};
    my %gene_pos_genome = %{$gene_pos{$genome}};
    my @sorted_genes = sort { $gene_pos_genome{$a} <=> $gene_pos_genome{$b} } keys %gene_pos_genome;
    foreach my $sample(@samples) {
	foreach my $gene(@sorted_genes) {
	    print OUTFILE $sdp,"\t",$sample,"\t",$gene,"\t",$gene_pos_genome{$gene},"\t",$gene_cov{$genome}{$gene}{$sample},"\n";
	}
    }
    close OUTFILE;
}
