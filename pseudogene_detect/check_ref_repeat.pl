#!/usr/bin/perl

######################################################################################################################################
## This program is to add the information of whether the candidate L1-mediated DEL is located at reference L1 region              ####
######################################################################################################################################

use strict;
use warnings;

my $rep_masker = shift; #hg19 repeat masker
my $query = shift; #somatic L1 candidates

die "perl $0 <repeat_masker> <somatic_ORFeus_Dox.list.txt>\n" unless ($rep_masker && $query);

my %h; my %h_line;
open RM, $rep_masker or die $!;
while(<RM>){
	chomp;
	my @a = split;
	my $chr = $a[5];
#	$chr =~ s/chr//;
	my $rep_type = $a[11];
	my $lpos = $a[6];
    my $rpos = $a[7];
	if($rep_type eq "LINE"){
		$h_line{$chr}{"$lpos-$rpos"} = "LINE";
	}elsif($rep_type eq "Simple_repeat"){
		$h{$chr}{"$lpos-$rpos"} = "Simple_repeat";
	}elsif($rep_type eq "Low_complexity"){
		$h{$chr}{"$lpos-$rpos"} = "Low_complexity";
	}elsif($rep_type eq "Satellite"){
		$h{$chr}{"$lpos-$rpos"} = "Other";
	}
}
close RM;

#print "L1_pos\tcns_lclip_cnt\tcns_rclip_cnt\tcnt_lpolyA_reads_from_sam_zx\tcnt_rpolyA_reads_from_sam_zx\tlpolyA_in_cns\trpolyA_in_cns\traw_ldisc_cnt\traw_rdisc_cnt\tcns_ldisc_cnt\tcns_rdisc_cnt\tref_repeat\tcnv_pos\tl_distance\tr_distance\tcnv_len\n";
open QU, $query or die $!;
while(<QU>){
	if(/pos/){
		print;
		next;
	}
    chomp;
    my @a = split;
	my $L1_pos = $a[1];
    my ($chr, $rclip_pos, $lclip_pos) = $L1_pos=~/(\w+):(\d+)-(\d+)/; #$rclip_pos < $lclip_pos

    my $ref_repeat = "nan";
    my $repeat_type = "nan";
	#for LINE	
	for my $pos_str(keys %{$h_line{$chr}}){
        my $rep_lpos = (split /-/, $pos_str)[0];
        my $rep_rpos = (split /-/, $pos_str)[1];
        if(abs($rclip_pos-$rep_lpos)<20 && abs($lclip_pos-$rep_rpos)<20){
            $repeat_type =  $h_line{$chr}{"$rep_lpos-$rep_rpos"};
            $ref_repeat = "$chr:$pos_str";
            last;
        }
	}
	#for other repeat elements
    for my $pos_str(keys %{$h{$chr}}){
		my $rep_lpos = (split /-/, $pos_str)[0];
		my $rep_rpos = (split /-/, $pos_str)[1];
#		print "$chr\t$rclip_pos\t$lclip_pos\t$rep_lpos\t$rep_rpos\n";
		# the 5' or 3' position of the repeat is located between the two endpoints of the repeat
		if(($rclip_pos>=$rep_lpos && $rclip_pos<=$rep_rpos) || ($lclip_pos>=$rep_lpos && $lclip_pos<=$rep_rpos)){
			$repeat_type =  $h{$chr}{"$rep_lpos-$rep_rpos"};
            $ref_repeat = "$chr:$pos_str";
            last;
		#the 5' or 3' position of the repeat is located close to either endpoints of the repeat
		}elsif(abs($rep_lpos-$rclip_pos)<10 || abs($rep_rpos-$lclip_pos)<10 || abs($rep_lpos-$lclip_pos)<10 || abs($rep_rpos-$rclip_pos)<10){
#		if(abs($lpos-$rep_lpos)<50 && abs($rpos-$rep_lpos)<50 ){
			$repeat_type =  $h{$chr}{"$rep_lpos-$rep_rpos"};
			$ref_repeat = "$chr:$pos_str";
			last;
		}
	}
	if(($ref_repeat eq "nan") || ($repeat_type eq "nan")){
#		print "$_\t$repeat_type;$ref_repeat\tpass\n";
		print "$_\n";
#	}else{
#		print "$_\t$repeat_type;$ref_repeat\tfiltered\n";
	}
}
close QU;
