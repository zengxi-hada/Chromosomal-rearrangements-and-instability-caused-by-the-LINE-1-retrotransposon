#!/usr/bin/perl

use strict;
use warnings;

my $ORFeus_disc_candidate = shift; #candidate_list_from_disc.txt, this file is actually not used in the program.
my $candidate_pre_flt = shift; #ORFeus_insertion.2-types.txt.add_clip_ORFeus_ratio.flt

die "perl $0 <candidate_list_from_disc.txt> <ORFeus_insertion.2-types.txt.add_clip_ORFeus_ratio.flt>\n" unless ($ORFeus_disc_candidate && $candidate_pre_flt);

#this block is no used.
my %h_disc;
open ODC, $ORFeus_disc_candidate or die $!;
while(<ODC>){
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $cns_lclip = $a[5]; #lclip reads map to ORFeus
	my $cns_rclip = $a[6]; #rclip reads map to ORFeus
	my $cns_ldisc = $a[8]; #ldisc reads map to ORFeus
	my $cns_rdisc = $a[9]; #rdisc reads map to ORFeus
	if($cns_ldisc == 0){
		$h_disc{right_disc}{"$chr\t$pos"} = "$cns_rclip\t$cns_rdisc";
	}elsif($cns_rdisc == 0){
		$h_disc{left_disc}{"$chr\t$pos"} = "$cns_lclip\t$cns_ldisc";
	}
}
close ODC;

#type1/2 chr     pos     lclip_cnt       rclip_cnt       cnt_mate_in_cns cns_lclip_cnt   cns_rclip_cnt   total_clip,ORFeus_clip(include sencondary align when also mapping to ORFeus)    cnt_polyA_reads_from_sam_zx     lpolyA_in_cns   rpolyA_in_cns   n_ldisc n_rdisc n_cns_ldisc     n_cns_rdisc
my %h_candidate;;
open CPF, $candidate_pre_flt or die $!;
while(<CPF>){
	next if /lclip_cnt/;
	chomp;
	my @a = split;
	my $type = $a[0];
	my $chr = $a[1];
	my $pos = $a[2];
	my $cns_lclip_cnt = $a[6];
	my $cns_rclip_cnt = $a[7];
	my $lpolyA_cnt = $a[10];
	my $rpolyA_cnt = $a[11];
	my $cns_ldisc_cnt = $a[14];
	my $cns_rdisc_cnt = $a[15];
	my $ratio_info = $a[8];
	my @ratio_arr = split /,/, $ratio_info;
	my $ratio = $ratio_arr[2];
	if($cns_rclip_cnt > 0){ #bk supported by both ORFeus disc reads and ORFeus clip reads
		$h_candidate{right_clip}{"$chr\t$pos"} = "$cns_rclip_cnt\t$rpolyA_cnt\t$cns_rdisc_cnt\t$ratio";
	}elsif($cns_lclip_cnt > 0){ #bk supported by both ORFeus disc reads and ORFeus clip reads
		$h_candidate{left_clip}{"$chr\t$pos"} = "$cns_lclip_cnt\t$lpolyA_cnt\t$cns_ldisc_cnt\t$ratio";
	}elsif($cns_rclip_cnt==0 && $cns_lclip_cnt==0){ #bk supported by only ORFeus disc reads
		$h_candidate{disc_only}{"$chr\t$pos"} = "$cns_lclip_cnt\t$lpolyA_cnt\t$cns_ldisc_cnt\t$ratio";
	}
}
close CPF; 

my %left_specific;
for my $lpos_str(keys %{$h_candidate{left_clip}}){
	my $flag = 0;
	my @l_arr = split /\t/, $lpos_str;
	my $lchr = $l_arr[0];
	my $lpos = $l_arr[1];
	my @l_cnt_arr = split /\t/, $h_candidate{left_clip}{$lpos_str};
	my $cns_lclip_cnt = $l_cnt_arr[0];
	my $lpolyA_cnt = $l_cnt_arr[1];
	my $cns_ldisc_cnt = $l_cnt_arr[2];
	my $l_ratio = $l_cnt_arr[3]; #ratio of clip reads mapping to cns
	for my $rpos_str(keys %{$h_candidate{right_clip}}){
		my @r_arr = split /\t/, $rpos_str;
		my $rchr = $r_arr[0];
        my $rpos = $r_arr[1];
		if(($lchr eq $rchr) && abs($lpos-$rpos)<=200){
			my @r_cnt_arr = split /\t/, $h_candidate{right_clip}{$rpos_str};
		    my $cns_rclip_cnt = $r_cnt_arr[0];
		    my $rpolyA_cnt = $r_cnt_arr[1];
		    my $cns_rdisc_cnt = $r_cnt_arr[2];
			my $r_ratio = $r_cnt_arr[3];
			my $type = $r_cnt_arr[4];
			print "$lchr:$lpos-$rpos\t$cns_lclip_cnt,$cns_rclip_cnt\t$lpolyA_cnt,$rpolyA_cnt\t$cns_ldisc_cnt,$cns_rdisc_cnt\t$l_ratio,$r_ratio\n"; #for both left and right
			$flag = 1;
		}
	}
#	if($flag==0 && $cns_lclip_cnt>=2){
	if($flag==0 && $cns_lclip_cnt>=1){ #updated at 2023-08-03
		$left_specific{$lpos_str} = $l_ratio;
	}
}

##for the insertions only have right or left side
#for insertions with both ORFeus clip and disc reads
my %right_specific;
for my $rpos_str(keys %{$h_candidate{right_clip}}){
	my $flag = 0;
	my @r_arr = split /\t/, $rpos_str;
	my $rchr = $r_arr[0];
	my $rpos = $r_arr[1];
	my @r_cnt_arr = split /\t/, $h_candidate{right_clip}{$rpos_str};
    my $cns_rclip_cnt = $r_cnt_arr[0];
    my $rpolyA_cnt = $r_cnt_arr[1];
    my $cns_rdisc_cnt = $r_cnt_arr[2];
	my $r_ratio = $r_cnt_arr[3];
    for my $lpos_str(keys %{$h_candidate{left_clip}}){
		my @l_arr = split /\t/, $lpos_str;
        my $lchr = $l_arr[0];
        my $lpos = $l_arr[1];
		if(($lchr eq $rchr) && abs($lpos-$rpos)<=200){
			$flag = 1;
		}
	}
#	if($flag==0 && $cns_rclip_cnt>=2){
	if($flag==0 && $cns_rclip_cnt>=1){ #updated at 2023-08-03
		$right_specific{$rpos_str} = $r_ratio;
	}
}

for my $pos_str(keys %left_specific){
	my @a = split /\t/, $pos_str;
	my $chr = $a[0];
	my $pos = $a[1];
	my $l_ratio = $left_specific{$pos_str};
	print "$chr:$pos-$pos\t$l_ratio\n";
}

for my $pos_str(keys %right_specific){
    my @a = split /\t/, $pos_str;
    my $chr = $a[0];
    my $pos = $a[1];
	my $r_ratio = $right_specific{$pos_str};
    print "$chr:$pos-$pos\t$r_ratio\n";
}

#for the insertions only have ORFeus reads
for my $pos_str(keys %{$h_candidate{disc_only}}){
	my $flag = 0;
    my @arr = split /\t/, $pos_str;
    my $chr = $arr[0];
    my $pos = $arr[1];
	print "$chr:$pos-$pos\tdisc_only\n";
}
