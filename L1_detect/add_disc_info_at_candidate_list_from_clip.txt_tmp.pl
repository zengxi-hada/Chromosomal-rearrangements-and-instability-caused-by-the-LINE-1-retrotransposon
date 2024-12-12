#!/usr/bin/perl

##################################################################################################################
## this program is to add number of left/right clipped reads in LINE-1 cns at candidate_list_from_disc.txt     ###
##################################################################################################################

use strict;
use warnings;

my $disc_sites = shift; #tmp/raw_discordant_reads_tmp0
my $clip_sites = shift; #candidate_list_from_clip.txt_tmp

die "perl $0 <tmp/raw_discordant_reads_tmp0> <candidate_list_from_clip.txt_tmp>\n" unless ($clip_sites && $disc_sites);

my %h;
open DS, $disc_sites or die $!;
while(<DS>){
    chomp;
    my @a = split;
    my $chr = $a[0];
    my $pos = $a[1];
	my $n_raw_ldisc = $a[2];
	my $n_raw_rdisc = $a[3];
	my $n_cns_ldisc = $a[4];
	my $n_cns_rdisc = $a[5];
	$h{"$chr-$pos"} = "$n_raw_ldisc\t$n_raw_rdisc\t$n_cns_ldisc\t$n_cns_rdisc";
}
close DS;

print "chr\tpos\tlclip_cnt\trclip_cnt\tcnt_mate_in_cns\tcns_lclip_cnt\tcns_rclip_cnt\tcnt_polyA_reads_from_sam_zx\tn_ldisc\tn_rdisc\tn_cns_ldisc\tn_cns_rdisc\n";
open CS, $clip_sites or die $!;
while(<CS>){
	next if /pos/;
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $lclip_cnt = $a[2];
    my $rclip_cnt = $a[3];
	my $cnt_mate_in_cns = $a[4];
	my $cns_lclip_cnt = $a[5];# reads whose left clipped part could map to LINE1-cns
    my $cns_rclip_cnt = $a[6];# reads whose right clipped part could map to LINE1-cns
	my $refined_polyA_reads = $a[7]; # get from sam, zengxi
	my ($cns_disc_info, $n_cns_ldisc, $n_cns_rdisc, $n_ldisc, $n_rdisc);
	if(exists $h{"$chr-$pos"}){
		$cns_disc_info = $h{"$chr-$pos"};
		$n_ldisc = (split /\t/, $cns_disc_info)[0]; #raw count
        $n_rdisc = (split /\t/, $cns_disc_info)[1]; # raw count
		$n_cns_ldisc = (split /\t/, $cns_disc_info)[2];
		$n_cns_rdisc = (split /\t/, $cns_disc_info)[3];
	}else{
		$n_ldisc = "merged_in_peak";
        $n_rdisc = "merged_in_peak";
		$n_cns_ldisc = "merged_in_peak";
		$n_cns_rdisc = "merged_in_peak";
	}
    print "$chr\t$pos\t$lclip_cnt\t$rclip_cnt\t$cnt_mate_in_cns\t$cns_lclip_cnt\t$cns_rclip_cnt\t$refined_polyA_reads\t$n_ldisc\t$n_rdisc\t$n_cns_ldisc\t$n_cns_rdisc\n";
}
close CS;
