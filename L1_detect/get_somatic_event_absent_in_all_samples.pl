#!/usr/bin/perl

use strict;
use warnings;

my $control_list = shift; #cns.txt.list, all sample
my $candidate_file = shift; #ORFeus_insertion.2-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic-5.2
my $somatic_events = shift; #ORFeus_insertion.2-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic_absent_in_all_samples-5.2

die "perl $0 <cns.txt.list> <ORFeus_insertion.2-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic-5.2> <ORFeus_insertion.2-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic-5.2_absent_from_all_samples>\n" unless ($control_list && $candidate_file && $somatic_events);

my @tmp_arr = split /\//, $candidate_file;
my $candidate_sample = $tmp_arr[0];

my %h_germline;
my %h_tmp;
open CTL, $control_list or die $!;
while(<CTL>){
	chomp;
	my $cns_txt = $_;
	my @tmp = split /\//, $cns_txt;
	my $sample = $tmp[0];
	next if $sample eq $candidate_sample;
	open CNS, $cns_txt or die $!;
	while(<CNS>){
	    next if /^#/;
	    chomp;
	    my @a = split;
	    my $chr = $a[0];
	    my $lpos = $a[2];
	    my $rpos = $a[3];
	    if($lpos == -1){
	        $lpos = $rpos;
	    }elsif($rpos == -1){
	        $rpos = $lpos;
	    }
		my @tmp_arr = ($lpos, $rpos);
        $rpos = &max(@tmp_arr);
        $lpos = &min(@tmp_arr);
#		$h_germline{$chr}{"$lpos\t$rpos"} = 1;
        $h_germline{$chr}{$lpos} = 1;
		$h_germline{$chr}{$rpos} = 1;
        my $mid = int(($lpos + $rpos) / 2);
        $h_tmp{$chr}{$mid} = 1;
	}
	close CNS;
}
close CTL;

open SME, ">$somatic_events" or die $!;
open CDF, $candidate_file or die $!;
while(<CDF>){
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $lpos = $a[1];
	my $rpos = $a[2];
	my $ori_lpos = $a[1];
	my $ori_rpos = $a[2];
	my @tmp_arr = ($lpos, $rpos);
    $rpos = &max(@tmp_arr);
    $lpos = &min(@tmp_arr);
	my $new_mid = int(($lpos + $rpos) / 2);
	my $flag = 0;
	for my $i(-25..25){
		my $new_pos = $new_mid + $i;
		if(exists $h_tmp{$chr}{$new_pos}){
			$flag = 1;
			last;
		}
	}
#	next if $flag == 1;

#	my ($chr, $lpos, $rpos) = $pos_str=~ /(\w+):(\d+)-(\d+)/;
#	if(not exists $h_germline{$chr}{"$lpos\t$rpos"}){ # a somatic event
	if($flag==0 || ((not exists $h_germline{$chr}{$lpos}) && (not exists $h_germline{$chr}{$rpos}))){ # a somatic event
		print SME "$chr\t$ori_lpos\t$ori_rpos\n";
	}
}
close CDF;
close SME;		

###############################subroutine##################
sub min
{
    my @a = @_;
    my $min = $a[0];
    for my $i(@a){
        next if $i == -1;
        if($i < $min){
            $min = $i;
        }
    }
    return $min;
}

sub max
{
    my @a = @_;
    my $max = $a[0];
    for my $i(@a){
        if($i > $max){
            $max = $i;
        }
    }
    return $max;
}
