#!/usr/bin/perl

use strict;
use warnings;

##############################################################################################################################################################
#this program is to calculate and add the ratio of reads mapping to ORFeus_L1 cns relative to the total number of reads covering on a specific genomic site  #
##############################################################################################################################################################

my $candidate = shift; #ORFeus_insertion.txt
my $sam = shift; #backup/*.bam.clipped.sam

die "perl $0 <ORFeus_insertion.txt>> <backup/*.bam.clipped.sam>" unless ($sam && $candidate);

my %h;
my %h_map_ratio;
open CAN, $candidate or die $!;
while(<CAN>){
	chomp;
	next if /cns_rclip_cnt/;
	my @a = split;
	my $chr = $a[1];
	my $pos = $a[2];
	$h{"$chr\t$pos"} = 1;
	$h_map_ratio{"$chr\t$pos"}{total_clip} = 0;
   	$h_map_ratio{"$chr\t$pos"}{ORFeus_clip} = 0;
}
close CAN;

open SAM, $sam or die $!;
while(<SAM>){
	next if /^@/;
	chomp;
	my @a = split;
	my $id = $a[0];
	my $cigar = $a[5];
	my @tmp = split /~/, $id;
	my $ori_chr = $tmp[3];
	my $ori_pos = $tmp[4];
	next if $cigar =~ /H/;
#	print "#$cigar\n";
	if(exists $h{"$ori_chr\t$ori_pos"}){
		my $seq = $a[9];
		my $len = length($seq);
		my $map_len = 0;
		while($cigar =~ /((\d+)[MSID])/g){
			my $align_len = $2;
			my $cigar_tmp = $1;
			if($1 =~ /M/){
				$map_len += $align_len;
#				print "#$cigar\t$cigar_tmp\t$align_len\n";
			}
		}
		my $map_ratio = sprintf("%.2f", $map_len/$len);
		$h_map_ratio{"$ori_chr\t$ori_pos"}{total_clip}++;
		if($map_ratio >= 0.6){
			$h_map_ratio{"$ori_chr\t$ori_pos"}{ORFeus_clip}++;
		}
	}
}
close SAM;

#chr     pos     lclip_cnt       rclip_cnt       cnt_mate_in_cns cns_lclip_cnt   cns_rclip_cnt   cnt_polyA_reads_from_sam_zx     lpolyA_in_cns   rpolyA_in_cns   n_ldisc n_rdisc n_cns_ldisc
#     n_cns_rdisc
open CAN, $candidate or die $!;
while(<CAN>){
	if(/cns_rclip_cnt/){
		print "type1/2\tchr\tpos\tlclip_cnt\trclip_cnt\tcnt_mate_in_cns\tcns_lclip_cnt\tcns_rclip_cnt\ttotal_clip,ORFeus_clip(include sencondary align when also mapping to ORFeus)\tcnt_polyA_reads_from_sam_zx\tlpolyA_in_cns\trpolyA_in_cns\tn_ldisc\tn_rdisc\tn_cns_ldisc\tn_cns_rdisc\n";
		next;
	}
    chomp;
    my @a = split;
    my $chr = $a[1];
    my $pos = $a[2];
	my $cns_lclip_cnt = $a[6];
	my $cns_rclip_cnt = $a[7];
	my $cns_clip_cnt = $cns_lclip_cnt + $cns_rclip_cnt;
	my $total_clip = $h_map_ratio{"$chr\t$pos"}{total_clip};
	my $ORFeus_clip = $h_map_ratio{"$chr\t$pos"}{ORFeus_clip};
	my $ratio;
	if($total_clip != 0){
		$ratio = sprintf("%.2f", $ORFeus_clip/$total_clip);
	}else{
		$ratio = 0;
	}
    print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$total_clip,$ORFeus_clip,$ratio\t$a[8]\t$a[9]\t$a[10]\t$a[11]\t$a[12]\t$a[13]\t$a[14]\n"
}
close CAN;
