#!/usr/bin/perl

use strict;
use warnings;

my $general_file = shift; #include cns_polyA reads number for each clip site; new candidate_list_from_clip.txt_tmp at L1_add_polyA_at_clip_stp
my $candidate_file = shift; #target file

die "perl $0 <general_file_include_polyA_info> <candidate_file>\n" unless ($general_file && $candidate_file); 

my %h;
open GF, $general_file or die $!;
while(<GF>){
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $lpolyA_in_cns = $a[7]; #number of left clip reads when maping to cns_polyA?
	my $rpolyA_in_cns = $a[8]; #number of right clip reads when maping to cns_polyA?
	my $polyA_in_cns = $lpolyA_in_cns + $rpolyA_in_cns;
	$h{"$chr-$pos"} = "$lpolyA_in_cns\t$rpolyA_in_cns";
}
close GF;

print "chr\tpos\tlclip_cnt\trclip_cnt\tcnt_mate_in_cns\tcns_lclip_cnt\tcns_rclip_cnt\tcnt_polyA_reads_from_sam_zx\tlpolyA_in_cns\trpolyA_in_cns\tn_ldisc\tn_rdisc\tn_cns_ldisc\tn_cns_rdisc\n";
open CF, $candidate_file or die $!;
while(<CF>){
	next if /pos/;
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $lclip_cnt = $a[2];
    my $rclip_cnt = $a[3];
	my $mate_in_cns = $a[4]; #number of clipped reads whose mates are mapped to LINE1-cns(consemsus)
    my $cns_lclip_cnt = $a[5]; # reads whose left clipped part could map to LINE1-cns
    my $cns_rclip_cnt = $a[6]; # reads whose right clipped part could map to LINE1-cns
	my $refined_polyA_reads = $a[7]; # get from sam, zengxi
	my $ldisc_cnt = $a[8];
	my $rdisc_cnt = $a[9];
    my $cns_ldisc_cnt = $a[10]; # reads whose left disc part could map to LINE1-cns
    my $cns_rdisc_cnt = $a[11]; # reads whose right disc part could map to LINE1-cns
	my $polyA_info = $h{"$chr-$pos"}; # polyA info obtained from xTea output: new candidate_list_from_clip.txt_tmp of call_TEI_candidate_sites_from_clip_reads_v2_mosaic()
	print "$chr\t$pos\t$lclip_cnt\t$rclip_cnt\t$mate_in_cns\t$cns_lclip_cnt\t$cns_rclip_cnt\t$refined_polyA_reads\t$polyA_info\t$ldisc_cnt\t$rdisc_cnt\t$cns_ldisc_cnt\t$cns_rdisc_cnt\n";
}
