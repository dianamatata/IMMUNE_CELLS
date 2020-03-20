#!/usr/bin/perl -w
use strict;
#Author: Cédric Howald
#November 15, 2017
#This script will create a filtered GTF file "output_prefix.gtf" and a metadata file listing which transcript_ids have been removed and why "output_prefix_rmTranscriptIds.tsv". Log information will be displayed on the standard output and should be redirected in a file (output_prefix.log)

#Usage: ./gencode_annotation_filtering.pl gencode.v19.annotation.gtf gencode.v19.annotation.filtered > gencode.v19.annotation.filtered.log

print STDOUT "This script will identify transcript ids of all non gene features (3rd field of GTF != \"gene\"), non protein coding, non lincRNA and non processed transcripts that overlap a protein coding or lincRNA gene (3rd field of GTF == \"gene\"). All non gene features having these transcript ids will be removed from the GTF file.\n\n";

my $bedtools = "/software/UHTS/Analysis/BEDTools/2.26.0/bin/intersectBed";

die "\nError: bedtools intersectBed is not executable: $bedtools\n\n" unless -X $bedtools;

my $input_GTF = $ARGV[0];
my $output_file = $ARGV[1];

my $out1=$output_file."_temp001.bed";
my $out2=$output_file."_temp002.bed";
open I, $input_GTF or die $!;
#parse the GTF and create a file containing the protein coding and lincRNA genes and another files containing non protein coding, lincRNA and processed transcripts features
open TMP1, ">$out1" or die $!;
open TMP2, ">$out2" or die $!;
print STDOUT "Parsing the GTF file...\n";
while(<I>){
    chomp(my $line = $_);
    if($line =~ /^#/){
        next;
    }
    my @l = split "\t", $line;
    my @info = split " ", $l[8];
    my $count = 0;
    my $transcript_id = "-";
    my $gene_type = "-";
    my $gene_id = "-";
    foreach(@info){
        my $tag = $_;
        if($tag eq "transcript_id"){
            $transcript_id = $info[$count+1];
            $transcript_id=~s/"//g;
            $transcript_id=~s/;//g;
        }elsif($tag eq "gene_type"){
            $gene_type = $info[$count+1];
            $gene_type=~s/"//g;
            $gene_type=~s/;//g;
        }elsif($tag eq "gene_id"){
            $gene_id = $info[$count+1];
            $gene_id=~s/"//g;
            $gene_id=~s/;//g;
        }
        $count++;
    }

    if($gene_type eq "-" || $gene_id eq "-"){
        `rm $out1 $out2`;
        die "\nError: the 9th field of the line doesn't contain the flag gene_type or gene_id.\nline: $line\n\n";
    }

    if($l[2] eq "gene"){
        if($gene_type eq "protein_coding" || $gene_type eq "lincRNA"){
            my $bed_start=$l[3]-1;
            print TMP1 $l[0]."\t". $bed_start."\t".$l[4]."\t".$gene_id."_".$gene_type."\n";
        }
    }else{
        if($gene_type ne "protein_coding" && $gene_type ne "processed_transcript" && $gene_type ne "lincRNA"){
            if($transcript_id eq "-"){
                `rm $out1 $out2`;
                die "\nError: the 9th field of the line doesn't contain the flag transcript_id.\nline: $line\n\n";
            }
            my $bed_start=$l[3]-1;
            print TMP2 $l[0]."\t". $bed_start."\t".$l[4]."\t".$l[2]."_".$gene_id."_".$transcript_id."_".$gene_type."\n";
        }
    }

}
close I;
close TMP1;
close TMP2;

#overlap the genes and features with bedTools
my $out3=$output_file."_temp003.bed";
`$bedtools -a $out1 -b $out2 -wa -wb > $out3`;

#extract the transcript ids that have to be removed from the original GTF file
print STDOUT "Extracting transcript ids that have to be removed...\n";
my %trans_id_to_remove;
open J, $out3 or die $!;
my $count = 0;
while(<J>){
    chomp(my $line = $_);
    my @tb = split "\t", $line;
    my @tb2 = split "_", $tb[7];
    my $trans_id = $tb2[2];
    my $startGene_1based=$tb[1]+1;
    my $startFeature_1based=$tb[5]+1;
    if(exists $trans_id_to_remove{$trans_id}){
        $trans_id_to_remove{$trans_id}.="/$tb[0]_$startGene_1based\_$tb[2]_$tb[3]<->$tb[4]_$startFeature_1based\_$tb[6]_$tb[7]";
    }else{
        $trans_id_to_remove{$trans_id}="$tb[0]_$startGene_1based\_$tb[2]_$tb[3]<->$tb[4]_$startFeature_1based\_$tb[6]_$tb[7]";
        $count++;
    }
}
close J;
print STDOUT "$count transcript ids have to be removed from the GTF...\n";

#create the output files containing the transcript id that were removed from the GTF with the overlap
my $output1 = $output_file.".rmTranscriptIds.tsv";
print STDOUT "Creating the file containing all transcript ids to removed ($output1)...\n";
open OUT1, ">$output1" or die $!;
print OUT1 "#This file contains all transcript ids removed from the GTF file (field 1) and why they have been removed (field 2). The field 2 show all overlaps (<->) between a gene and a feature having the transcript id of field 1 (gene1_info<->feature1_info/gene2_info<->feature2_info/gene3_info<->feature3_info/...)\n";
foreach my $k (keys %trans_id_to_remove){
    print OUT1 $k."\t".$trans_id_to_remove{$k}."\n"
}
close OUT1;

#create the filtered GTF file
my $output2 = $output_file.".gtf";
print STDOUT "Creating the filtered GTF file ($output2)...\n";
my %feature_removed;
open I, $input_GTF or die $!;
open OUT2, ">$output2" or die $!;
print OUT2 "##Applied filter: this GTF file has been filtered for non gene features (3rd field) that are non protein_coding, non processed_transcript and non lincRNA (gene_type field) overlapping a gene feature (3rd field = \"gene\") that is protein_coding or lincRNA. All features except genes matching a transcript id detected as having a such overlap will be removed from the GTF file\n";

while(<I>){
    chomp(my $line = $_);
    if($line =~ /^#/){
        print OUT2 $line."\n";
        next;
    }
    my @l = split "\t", $line;
    my @info = split " ", $l[8];
    my $count = 0;
    my $transcript_id = "-";
    foreach(@info){
        my $tag = $_;
        if($tag eq "transcript_id"){
            $transcript_id = $info[$count+1];
            $transcript_id=~s/"//g;
            $transcript_id=~s/;//g;
        }
        $count++;
    }

    if($l[2] eq "gene"){
        print OUT2 $line."\n";
    }else{
        if($transcript_id eq "-"){
            `rm $out1 $out2 $out3 $output1 $output2`;
            die "\nError: the 9th field of the line doesn't contain the flag transcript_id.\nline: $line\n\n";
        }
        if(!exists $trans_id_to_remove{$transcript_id}){
            print OUT2 $line."\n";
        }else{
            $feature_removed{$l[2]}++;            
        }
    }
}
close I;
close OUT2;

print STDOUT "\nfeatures filtered out (3rd field of the GTF). The total should correspond to the number of lines removed:\n";
for my $k (sort keys %feature_removed){
    print "-$k\t$feature_removed{$k}\n";
}

`rm $out1 $out2 $out3`;
print STDOUT "\nFiltering done!\n\n";

