#! /usr/bin/perl -w
use Getopt::Long;

#Usage: perl core_cov.pl orthofile 

my @samples = qw(F01 F02 F03 F04 F06 F07 F08 F09 F10 F11 F12 F13 F14 F15 F16 N01 N02 N03 N04 N06 N07 N08 N09 N10 N11 N12 N13 N14 N15 N16);

my %SDP_members = (
'C4S75' => 'api_1',
'Ga0307799' => 'api_1',
'Ga0312307' => 'api_1',
'C4S76' => 'api_apis_dorsa',
'C4S77' => 'api_apis_dorsa',
'Ga0061079' => 'api_bombus',
'BBC0122' => 'bapis', 
'BBC0178' => 'bapis', 
'BBC0244' => 'bapis',
'H3S78' => 'bapis',
'H3S79' => 'bapis',
'H3S80' => 'bapis',
'H3S81' => 'bapis',
'H3S82' => 'bapis',
'H3S83' => 'bapis',
'H3S84' => 'bapis',
'H3S85' => 'bapis',
'H3S86' => 'bapis',
'H3S89' => 'bapis',
'H3S90' => 'bapis',
'H3S94' => 'bapis',
'H3S95' => 'bapis',
'H3S98' => 'bapis',
'H3U94' => 'bapis',
'H3V12' => 'bapis',
'H3V13' => 'bapis', 
'H3V14' => 'bapis', 
'H3V15' => 'bapis', 
'H3V16' => 'bapis', 
'H3V17' => 'bapis', 
'PEB0122' => 'bapis',
'PEB0149' => 'bapis',
'PEB0150' => 'bapis',
'BAST' => 'bifido_1.3',
'Ga0072395' => 'bifido_1.1',
'Ga0072396' => 'bifido_1.2',
'Ga0072401' => 'bifido_1.4',
'Ga0133549' => 'bifido_1.3',
'Ga0133551' => 'bifido_1.4',
'Ga0133552' => 'bifido_1.4',
'Ga0133553' => 'bifido_1.5',
'Ga0312305' => 'bifido_1_cerana',
'H3S92' => 'bifido_1.3',
'H3S93' => 'bifido_1.3',  
'H3S97' => 'bifido_1.2', 
'H3T81' => 'bifido_1.1',
'H3T86' => 'bifido_1.2', 
'H3T89' => 'bifido_1.3', 
'H3T90' => 'bifido_1.4', 
'H3T91' => 'bifido_1.4', 
'H3T92' => 'bifido_1.1', 
'H3T93' => 'bifido_1.4', 
'H3U87' => 'bifido_1.2', 
'H3U88' => 'bifido_1.1', 
'H3U89' => 'bifido_1.2', 
'H3U98' => 'bifido_1.2',
'H3V04' => 'bifido_1.2', 
'H3V06' => 'bifido_1.4',  
'BCOR' => 'bifido_2',
'BINDI' => 'bifido_2',
'Ga0072398' => 'bifido_2',
'Ga0133550' => 'bifido_2', 
'H3S87' => 'bifido_2', 
'H3V08' => 'bifido_2', 
'Ga0056940' => 'bifido_bombus',
'Ga0057557' => 'bifido_bombus',
'Ga0098206' => 'bifido_bombus',
'Ga0372762' => 'bom_1', 
'Ga0216355' => 'bom_1', 
'Ga0216358' => 'bom_1', 
'Ga0216357' => 'bom_apis_melli', 
'Ga0308518' => 'bom_bombus',
'D9V35' => 'com_1', 
'Ga0216351' => 'com_1',
'Ga0216349' => 'com_1',
'Ga0216352' => 'com_1',
'Ga0216354' => 'com_1',
'Ga0216356' => 'com_1',
'Ga0216359' => 'com_1',
'Ga0216360' => 'com_1',
'H3R28' => 'com_1',
'H3R30' => 'com_1',
'H3R31' => 'com_1',
'H3R32' => 'com_1',
'H3R33' => 'com_1',
'H3R39' => 'com_1',
'H3S70' => 'com_1',
'H3S71' => 'com_1',
'H3S72' => 'com_1',
'H3T44' => 'com_1',
'H3T46' => 'com_1',
'H3T50' => 'com_1',
'H3T52' => 'com_1',
'H3T57' => 'com_1',
'H3T58' => 'com_1',
'H3T59' => 'com_1',
'H3V05' => 'com_1',
'CIN' => 'com_drosophila', 
'Ga0248239' => 'com_monarch',
'Ga0070888' => 'firm4_1',
'Ga0227356' => 'firm4_1',
'Ga0227357' => 'firm4_1',
'Ga0227358' => 'firm4_1',
'Ga0227360' => 'firm4_1',
'Ga0326456' => 'firm4_1',
'H3T38' => 'firm4_1',
'H3T39' => 'firm4_1', 
'H3T40' => 'firm4_1', 
'H3T41' => 'firm4_1', 
'H3T42' => 'firm4_1', 
'H3T51' => 'firm4_1', 
'Ga0072399' => 'firm4_2',
'H3R26' => 'firm4_2',
'Ga0226842' => 'firm5_1',
'Ga0072402' => 'firm5_1',
'Ga0061073' => 'firm5_1',
'Ga0133563' => 'firm5_1',
'H3R21' => 'firm5_1',
'H3R25' => 'firm5_1',
'H3T48' => 'firm5_1',
'Ga0133561' => 'firm5_2',
'Ga0072400' => 'firm5_2',
'LACWKB8' => 'firm5_2',
'H3T47' => 'firm5_2', 
'Ga0226840' => 'firm5_3',
'Ga0227359' => 'firm5_3',
'Ga0133562' => 'firm5_3',
'Ga0072404' => 'firm5_3',
'Ga0225908' => 'firm5_3',
'H3R22' => 'firm5_3',
'H3R27' => 'firm5_3',
'H3R29' => 'firm5_3',
'Ga0070887' => 'firm5_4',
'Ga0226839' => 'firm5_4',
'Ga0133564' => 'firm5_4',
'LACWKB10' => 'firm5_4',
'Ga0072403' => 'firm5_4',
'H3T43' => 'firm5_4',
'H3U40' => 'firm5_4',
'H3U49' => 'firm5_4',
'H3U50' => 'firm5_4',
'Ga0226852' => 'firm5_bombus',
'Ga0226847' => 'firm5_bombus',
'Ga0226843' => 'firm5_bombus',
'Ga0307804' => 'firm5_7',
'Ga0312303' => 'firm5_7',
'Ga0311394' => 'firm5_7',
'Ga0312304' => 'firm5_7',
'Ga0133554' => 'fper',
'Ga0077910' => 'fper',
'A9G14' => 'gilli_1',
'A9G17' => 'gilli_1', 
'B6D03' => 'gilli_1',
'B6D08' => 'gilli_1',
'B6D19' => 'gilli_1',
'B6D22' => 'gilli_1',
'B5798' => 'gilli_1', 
'B5799' => 'gilli_1',
'B5800' => 'gilli_1',
'B5801' => 'gilli_1', 
'B5802' => 'gilli_1', 
'B5803' => 'gilli_1', 
'B5804' => 'gilli_1', 
'B5S40' => 'gilli_1', 
'B5S41' => 'gilli_1',
'B5S42' => 'gilli_1',
'B5S43' => 'gilli_1',
'B5S44' => 'gilli_1',
'B6C87' => 'gilli_1',
'B6C91' => 'gilli_1', 
'B6D02' => 'gilli_1', 
'B6D04' => 'gilli_1', 
'B6D05' => 'gilli_1', 
'B6D07' => 'gilli_1', 
'B6D11' => 'gilli_1', 
'B6D12' => 'gilli_1', 
'B6D13' => 'gilli_1', 
'B6D15' => 'gilli_1', 
'B6D16' => 'gilli_1', 
'B6D20' => 'gilli_1', 
'Ga0133559' => 'gilli_1',
'GAPWK' => 'gilli_1', 
'H3S73' => 'gilli_1', 
'H3S75' => 'gilli_1',
'H3S76' => 'gilli_1',
'H3S88' => 'gilli_1',
'H3S91' => 'gilli_1',
'H3T96' => 'gilli_1',
'A9G15' => 'gilli_2',
'A9G16' => 'gilli_2',
'A9G19' => 'gilli_2',
'B6C84' => 'gilli_2',
'B6C86' => 'gilli_2', 
'B6C88' => 'gilli_2',
'B6C89' => 'gilli_2',
'B6C90' => 'gilli_2',
'B6C92' => 'gilli_2', 
'B6C94' => 'gilli_2', 
'B6C96' => 'gilli_2',
'B6C97' => 'gilli_2', 
'B6C98' => 'gilli_2', 
'B6C99' => 'gilli_2',
'B6D06' => 'gilli_2', 
'B6D09' => 'gilli_2', 
'B6D10' => 'gilli_2', 
'B6D14' => 'gilli_2',
'B6D17' => 'gilli_2',
'B6D18' => 'gilli_2',
'B6D21' => 'gilli_2',
'B6D23' => 'gilli_2', 
'B6D26' => 'gilli_2', 
'Ga0133555' => 'gilli_2',
'Ga0133556' => 'gilli_2',
'Ga0227303' => 'gilli_2', 
'H3S74' => 'gilli_2', 
'H3S77' => 'gilli_2',
'H3T54' => 'gilli_2',
'H3T56' => 'gilli_2',
'H3T60' => 'gilli_2',
'H3T61' => 'gilli_2',
'H3T80' => 'gilli_2',
'H3T83' => 'gilli_2',
'H3T85' => 'gilli_2',
'H3U96' => 'gilli_2',
'H3U97' => 'gilli_2',
'Ga0133557' => 'gilli_3',
'Ga0133558' => 'gilli_3', 
'A9G10' => 'gilli_4',
'Ga0307802' => 'gilli_4',
'Ga0307800' => 'gilli_5',
'Ga0307801' => 'gilli_5', 
'A9G07' => 'gilli_6',
'A9G08' => 'gilli_6',
'Ga0307803' => 'gilli_6', 
'A9G03' => 'gilli_apis_andre',
'A9G11' => 'gilli_apis_dorsa',
'A9G12' => 'gilli_apis_dorsa',
'A9G13' => 'gilli_apis_dorsa',
'A9G23' => 'gilli_bombus',
'A9G24' => 'gilli_bombus',
'A9G29' => 'gilli_bombus',
'A9G30' => 'gilli_bombus',
'A9G32' => 'gilli_bombus',
'A9G35' => 'gilli_bombus',
'A9G37' => 'gilli_bombus',
'A9G38' => 'gilli_bombus',
'A9G39' => 'gilli_bombus',
'A9G41' => 'gilli_bombus',
'A9G42' => 'gilli_bombus',
'A9G43' => 'gilli_bombus',
'A9G44' => 'gilli_bombus',
'A9G47' => 'gilli_bombus',
'Ga0055041' => 'gilli_bombus',
'Ga0061080' => 'gilli_bombus',
'Ga0227298' => 'gilli_bombus',
'Ga0227301' => 'gilli_bombus',
'Ga0227302' => 'gilli_bombus',
'BGH94' => 'snod_1',
'BGH95' => 'snod_1',
'BGH96' => 'snod_1',
'BGH97' => 'snod_1',
'BGH98' => 'snod_1',
'BGH99' => 'snod_1',
'BGI00' => 'snod_1',
'BGI01' => 'snod_1',
'BGI02' => 'snod_1',
'BGI03' => 'snod_1',
'BGI04' => 'snod_1',
'BGI05' => 'snod_1',
'BGI06' => 'snod_1',
'BGI07' => 'snod_1',
'BGI08' => 'snod_1',
'BGI09' => 'snod_1',
'BGI10' => 'snod_1',
'BGI11' => 'snod_1',
'BGI12' => 'snod_1',
'BGI13' => 'snod_1',
'BGI14' => 'snod_1',
'BGI15' => 'snod_1',
'BHC42' => 'snod_1',
'BHC45' => 'snod_1',
'BHC47' => 'snod_1',
'BHC50' => 'snod_1',
'BHC52' => 'snod_1',
'BHC56' => 'snod_1',
'CPT77' => 'snod_1', 
'Ga0227304' => 'snod_1',
'Ga0227305' => 'snod_1',
'Ga0227306' => 'snod_1',
'Ga0326500' => 'snod_1',
'H3T75' => 'snod_1',
'H3T77' => 'snod_1',
'H3T79' => 'snod_1',
'H3T82' => 'snod_1',
'H3T84' => 'snod_1',
'H3U73' => 'snod_1',
'H3U82' => 'snod_1',
'H3V00' => 'snod_1', 
'H3V01' => 'snod_1', 
'H3V10' => 'snod_1', 
'H3V11' => 'snod_1', 
'SALWKB2' => 'snod_1',  
'BHC43' => 'snod_2',
'BHC53' => 'snod_2',
'BHC54' => 'snod_2',
'BGI32' => 'snod_bombus',
'BGI34' => 'snod_bombus',
'BGI36' => 'snod_bombus',
'BHC44' => 'snod_bombus',
'BHC46' => 'snod_bombus',
'BHC48' => 'snod_bombus',
'BHC57' => 'snod_bombus',
'Ga0227310' => 'snod_bombus',
'SALWKB12' => 'snod_bombus',
'Ga0098617' => 'lkun',
'Ga0308523' => 'lkun',
'H3T37' => 'lkun',
);

my %ref_ids = (  
'firm5_1' => 'Ga0133563',
'firm5_2' => 'Ga0133561',
'firm5_3' => 'Ga0133562',
'firm5_4' => 'Ga0133564',
'firm5_bombus' => 'Ga0226852',
'firm5_6' => 'Ga0226847',
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
'bifido_1_cerana' => 'Ga0312305', 
'bifido_2' => 'BINDI',
'bifido_bombus' => 'Ga0098206',
'firm4_1' => 'Ga0326456', 
'firm4_2' => 'Ga0072399',
'fper' => 'Ga0077910',
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

#Read ortho-file, save in hash. Get all genome ids. 

my %OG_fam;
my %genome_ids;
my %gene_OG;
open(FILE, $ARGV[0]) or
    die "Cant open $ARGV[0]!";
while(<FILE>) {
    chomp;
    my @split = split(" ", $_);
    my $OG = shift @split;
    my $OG_id = substr $OG, 0, -1;
    foreach my $gene(@split) {
	$gene_OG{$gene} = $OG_id;
	my @gene_split  = split("_", $gene);
	my $genome_id = $gene_split[0];
	my $genome_sdp = $SDP_members{$genome_id};
	$genome_ids{$genome_sdp}{$genome_id}=1;
	$OG_fam{$genome_sdp}{$OG_id} = [] unless exists $OG_fam{$genome_sdp}{$OG_id};
	push @{$OG_fam{$genome_sdp}{$OG_id}},$gene;
    }
}
close FILE;

#Read the bed-files corresponding to ref_ids, and get the reference genome positions for each OG.
my %OG_ref_pos;
foreach my $sdp(keys %OG_fam) {
    my $ref_bed = $ref_ids{$sdp}.".bed";
    open(FILE, $ref_bed) or
	die "Cant open $ref_bed!";
    while(<FILE>) {
	chomp;
	my @split = split("\t",$_);
	my $gene = $split[3];
	if (exists $gene_OG{$gene}) {
	    my $genefam = $gene_OG{$gene};
	    my $start_pos = $split[1];	
	    $OG_ref_pos{$sdp}{$genefam} = $start_pos;
	}
    }
}

#Get the coverage of all gene-families in orthofile, for all listed samples
my %genefam_cov;
foreach my $sample(@samples) {
    my $bam_file = $sample."_vs_db_filt2x_sorted.bam";
    print "Processing bam-file: ",$bam_file,"\n";
    foreach my $sdp( sort keys %genome_ids) {
	my %genomes = %{$genome_ids{$sdp}};
	my @genomes = sort keys %genomes;
	my %gene_cov_sdp=();
	foreach my $genome(@genomes) {
	    my $bed_file = $genome.".bed";
	    my @cov_lines = `samtools bedcov $bed_file $bam_file`;
	    foreach my $line(@cov_lines) {
		chomp($line);
		my @split_line = split("\t",$line);
		my $gene_id = $split_line[3];
		my $gene_length = $split_line[2]-$split_line[1];
		my $gene_cov = sprintf("%.2f",($split_line[4]/$gene_length));
		$gene_cov_sdp{$gene_id} = $gene_cov;
	    }
	}
	my %OG_fam_sdp = %{$OG_fam{$sdp}};
	foreach my $OG(keys %OG_fam_sdp) {
	    my @genes = @{$OG_fam_sdp{$OG}};
	    my $genefam_cov = 0;
	    foreach my $gene(@genes) {
		$genefam_cov += $gene_cov_sdp{$gene};
	    }
	    $genefam_cov{$sample}{$sdp}{$OG} = $genefam_cov;
	}
    }
}

#Generate output file.


my @split_name = split("_",$ARGV[0]);
my $outfile = $split_name[0]."_corecov.txt";
open OUTFILE, '>'.$outfile or die $!;

print OUTFILE "SDP\tSample\tOG\tRef_pos\tCoverage\n";
foreach my $sdp( sort keys %genome_ids) {
    my %sdp_ref_pos = %{$OG_ref_pos{$sdp}};
    my @sorted_OG = sort { $sdp_ref_pos{$a} <=> $sdp_ref_pos{$b} } keys %sdp_ref_pos;
    foreach my $sample(@samples) {
	foreach my $OG(@sorted_OG) {
	    print OUTFILE $sdp,"\t",$sample,"\t",$OG,"\t",$sdp_ref_pos{$OG},"\t", $genefam_cov{$sample}{$sdp}{$OG},"\n";;
	}
    }
}
