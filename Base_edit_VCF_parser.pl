#!/usr/bin/perl -w

use strict;

# Data need to be analysed as follows
# 1. Trim adaptors (e.g. trim_galore)
# 2. Align with Bowtie2:
#       bowtie2 -p 4 -X 1000 -x /databank/indices/bowtie2/hg19/hg19 -1 filename_R1.fq -2 filename_R2.fq  -S filename.sam
# 3. Use samtools to call sequences
#       samtools view -S -b -o filename.bam filename.sam
#       samtools sort filename.bam  -o $filename.sorted.bam
#       samtools mpileup -f hg19.fa -u -v -o filename.vcf filename.sorted.bam
# 4. Use vcftools to extract coordinates of editing windows using bed files of the positive and negative stranded off-target editing windows
#       vcftools --vcf filename.vcf --bed Positive_strand_coords.bed --out ${filename}_Pos_strand.txt --recode --recode-INFO-all
#       vcftools --vcf filename.vcf --bed Negative_strand_coords.bed --out ${filename}_Neg_strand.txt --recode --recode-INFO-all
# 5. Analyse the resulting files with this script

my @files = $ARGV[0]

my %counters;
my %allele;
my %allele_freq;
my %list;
open(FHOUT, ">VCFparser_multi_ouptut_$strand.txt") or die "Cannot open file output file!\n";
open(FHOUTLIST, ">VCFparser_list_$strand.txt") or die "Cannot open file output file!\n";

foreach my $filename_and_path(@files)
{
unless ($filename_and_path =~ /(.*).txt.recode.vcf/) {die"filename does not match vcf format"};
my $file_name=$1;
my $file_path="";
if ($file_name =~ /(.*\/)(\V++)/) {$file_path = $1; $file_name = $2};

my $filename_out= $file_path.$file_name."_".$strand."_allele_counts.txt";

open(VCFFH, $filename_and_path) or die "Cannot open file $filename_and_path $!\n";
my $headder = "reference\tfirst_alternate\tdepth\talternate_num\tfrequency\tbase_edit_top\tbase_edit_bottom\tthreshold_base_edit_top\tthreshold_base_edit_bottom\tlikely_het";

print FHOUTLIST "file\t$headder\n";

while (my $line = <VCFFH>)
  {
    chomp $line;
    my @data = split(/\t/, $line);
    if ($line =~ /^chr/)
    {
        $counters{"lines_read"}++;
        #print "\n".$line;
        my $chr = $data[0]; if ($chr =~ /chr(.*)/){$chr = $1};
        my $position = $data[1];
        my $reference=$data[3];
        my $alternate=$data[4];
        my $first_alternate=substr($data[4],0,1);
        my $info = $data[7]; 
        my $depth = 0;
        my $allele_freq = 0;
       
        if($line =~ /.*DP=(\d++).*/){$depth = $1; $counters{"depth_detemined"}++;}
        else {print "$info\n"}
        if (($depth > 1000))
            {
            my $I_string;
            my $ref_for;
            my $ref_rev;
            my $non_for;
            my $non_rev;
            if($info =~ /.*I16=(.*);(.*)/)
                {
                $I_string = $1;
                ($ref_for, $ref_rev, $non_for, $non_rev, my $rest) = split (/,/, $I_string);
                my $denominator = $ref_for + $ref_rev + $non_for + $non_rev;
                my $alternate_num=$non_for+$non_rev;
                if ($denominator != 0){$allele_freq = ($ref_for + $ref_rev)/($denominator)};
                my $frequency = ($non_for+$non_rev)/$depth;
                $allele{$reference.$first_alternate}{"00count"}++;
                $allele{$reference.$first_alternate}{"01depth"}+=$depth;
                $allele{$reference.$first_alternate}{"02reference"}+=$ref_for+$ref_rev;
                $allele{$reference.$first_alternate}{"03alternate"}+=$non_for+$non_rev;
                
                my $base_edit_top=0; if (($reference eq "A") and ($first_alternate="G")){$base_edit_top=1}
                my $base_edit_bottom=0; if (($reference eq "T") and ($first_alternate="C")){$base_edit_bottom=1}

                my $threshold_edit=0; my $threshold_base_edit_top=0; my $threshold_base_edit_bottom=0;
                if (($frequency > 0.001) and ($frequency <0.1))
                {
                    $allele{$reference.$first_alternate}{"04_mutant_count"}++;
                    if ($base_edit_top==1){$threshold_base_edit_top=1}
                    if ($base_edit_bottom==1){$threshold_base_edit_bottom=1}
                };
                
                
                my$likely_het =0;
                if ($frequency > 0.0999){$allele{$reference.$first_alternate}{"05_het_count"}++};

                my $string = "$data[0]\t$data[1]";
                $list{$string}{$file_name}="$reference\t$first_alternate\t$depth\t$alternate_num\t$frequency\t$base_edit_top\t$base_edit_bottom\t$threshold_base_edit_top\t$threshold_base_edit_bottom\t$likely_het";
                print FHOUTLIST "$file_name\t$string\t$reference\t$first_alternate\t$depth\t$alternate_num\t$frequency\t$base_edit_top\t$base_edit_bottom\t$threshold_base_edit_top\t$threshold_base_edit_bottom\t$likely_het\n"
                }
            }
    }
    $counters{"Total_lines"}++;
  }
}

output_hash(\%counters);
output_2Dhash(\%allele, \*FHOUT);
print FHOUT "\n\nSite table\n";
output_2Dhash_special(\%list, \*FHOUT);

#This ouputs a 2column hash to a file
sub output_hash
{
    my ($hashref) = @_;
    foreach my $value (sort keys %$hashref)
    {
    print "$value\t".$$hashref{$value}."\n";
    }        
}

#This ouputs a 2D hash in a tab separated spreadsheet format
sub output_2Dhash
{
    my ($hashref, $filehandleout_ref) = @_;
        my %columns;
    foreach my $row (sort {$a cmp $b} keys %$hashref)
    {
                foreach my $column (sort {$a cmp $b} keys %{$$hashref{$row}}){$columns{$column}++} # gets a list of all of the columns in the file
    }
        my @columns = sort {$a cmp $b} keys %columns;
        foreach my $column (@columns){print $filehandleout_ref "\t$column"}
        foreach my $row (sort {$a cmp $b} keys %$hashref)
    {
                print $filehandleout_ref "\n$row";
                foreach my $column (@columns)
                {
                        if (exists$$hashref{$row}{$column}) 
            {
                if ($$hashref{$row}{$column} ne "")
                {
                print $filehandleout_ref "\t".$$hashref{$row}{$column}
                }
                else {print $filehandleout_ref "\t0"}
            }
                        else {print $filehandleout_ref "\t0"}
                }
        }
}

#This ouputs a 2D hash in a tab separated spreadsheet format
sub output_2Dhash_special
{
    my ($hashref, $filehandleout_ref) = @_;
    my $headder = "reference\tfirst_alternate\tdepth\talternate_num\tfrequency\tbase_edit_top\tbase_edit_bottom\tthreshold_base_edit_top\tthreshold_base_edit_bottom\tlikely_het";
        my %columns;
    foreach my $row (sort {$a cmp $b} keys %$hashref)
    {
                foreach my $column (sort {$a cmp $b} keys %{$$hashref{$row}}){$columns{$column}++} # gets a list of all of the columns in the file
    }
    #print Dumper (\%columns);
        my @columns = sort {$a cmp $b} keys %columns;
    #print Dumper (\@columns);
    print $filehandleout_ref "chr\tposition";
        foreach my $column (@columns){print $filehandleout_ref "\t$column$headder"}
        foreach my $row (sort {$a cmp $b} keys %$hashref)
    {
                print $filehandleout_ref "\n$row";
                foreach my $column (@columns)
                {
                        if (exists$$hashref{$row}{$column}) 
            {
                if ($$hashref{$row}{$column} ne "")
                {
                print $filehandleout_ref "\t".$$hashref{$row}{$column}
                }
                else {print $filehandleout_ref "\t0"}
            }
                        else {print $filehandleout_ref "\t-\t-\t0\t0\t0\t0\t0\t0\t0\t0"}
                }
        }
}
