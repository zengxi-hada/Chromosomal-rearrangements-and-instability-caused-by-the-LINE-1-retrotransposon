#!/usr/bin/perl

use strict;
use warnings;

####################################################################################################################################################
## This program is to select the high-confident pseudogene candiates by requiring the clip parts / disc mates of both sides map to close locations #
####################################################################################################################################################

my $candidate_pseudogenes = shift; #pseudogenes_rm_rep.txt_flt_by_near_clip

die "perl $0 <pseudogenes_rm_rep.txt_flt_by_near_clip>\n" unless $candidate_pseudogenes;

open CP, $candidate_pseudogenes or die $!;
while(<CP>){
	if(/clip_cnt/ || /only_left/ || /only_right/){
		print;
		next;
	}
	chomp;
	my @a = split;
	my $pos_str = $a[1];
	my $two_most_frq_supple = $a[5];
	my $clip_mates = $a[7];
	my $disc_mates = $a[8];
	my @frq_supple = split /;/, $two_most_frq_supple;
	my @left_frq_supple = split /,/, $frq_supple[0];
	my @right_frq_supple = split /,/, $frq_supple[1];
	my @disc_mates_arr = split /;/, $disc_mates;
	my @left_disc_mates_arr = split /,/, $disc_mates_arr[0];
	my @right_disc_mates_arr = split /,/, $disc_mates_arr[1];
	my ($flag_supple, $flag_disc_mate) = (0, 0);
	for my $i1(@left_frq_supple){
		my ($chr1, $pos1) = split /:/, $i1;
		for my $i2(@right_frq_supple){
			my ($chr2, $pos2) = split /:/, $i2;
			if(($chr1 eq $chr2) && abs($pos1-$pos2)<5000){
				$flag_supple = 1;
			}
		}
	}

	for my $i1(@left_disc_mates_arr){
        my ($chr1, $pos1) = split /:/, $i1;
        for my $i2(@right_disc_mates_arr){
            my ($chr2, $pos2) = split /:/, $i2;
            if(($chr1 eq $chr2) && abs($pos1-$pos2)<10000){
                $flag_disc_mate = 1;
            }
        }
    }
	#select the candidates whose supplementary align of left and right clip reads are close, and the mate of left / right disc reads map close
	if($flag_supple==1 && $flag_disc_mate==1){
		print "$_\n";
	}
}
close CP;
