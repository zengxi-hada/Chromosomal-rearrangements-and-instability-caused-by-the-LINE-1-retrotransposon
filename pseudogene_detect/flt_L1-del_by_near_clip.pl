#!/usr/bin/perl

use strict;
use warnings;

###############################################################################################################################################
## If there are another clip site with the same direction as the tranlocation candidate site in nearby region (±50bps),                     
## then this is considered at a tranlocation, which will be filtered out.      
###########################################################################################################################################

my $L1bk = shift; #transfered candidate_list_from_clip.txt_tmp
my $candidate_pseudogene = shift; #pseudogene.txt

die "perl $0 <candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc.add_xTea_polyA> <pseudogene.txt>\n" unless ($L1bk && $candidate_pseudogene);

my %h_bk;
open L1BK, $L1bk or die $!; #from xTea
while(<L1BK>){
	next if /pos/;
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $lclip_cnt = $a[2];
	my $rclip_cnt = $a[3];
	my $cns_lclip_cnt = $a[5]; # reads whose left clipped part could map to LINE1-cns
	my $cns_rclip_cnt = $a[6]; # reads whose right clipped part could map to LINE1-cns
	my $polyA_from_sam = $a[7]; # polyA reads from sam, zengxi
	my $lpolyA_in_cns = $a[8]; # number of left clip reads when maping to cns_polyA? from xTea
	my $rpolyA_in_cns = $a[9]; # number of right clip reads when maping to cns_polyA?
	my $ldisc_cnt = $a[10]; # raw count
    my $rdisc_cnt = $a[11]; # raw count
	my $cns_ldisc_cnt = $a[12]; # reads whose left disc part could map to LINE1-cns
    my $cns_rdisc_cnt = $a[13]; # reads whose right disc part could map to LINE1-cns
#	next if ($cns_ldisc_cnt eq "merged_in_peak" || $cns_rdisc_cnt eq "merged_in_peak"); # it means this site is merge into peaks beacuse it is close to other sites (merged peak: tmp/clip_peak_candidate.list)
	if($lclip_cnt  == 0){
#		if($cns_rclip_cnt>0 || $cns_ldisc_cnt>0 || $cns_rdisc_cnt>0){  #right part of a read was clipped (left side of a LINE-1)
			$h_bk{right_part_clipped}{$chr}{$pos} = $rclip_cnt;
#		}
	}elsif($rclip_cnt == 0){
#		if($cns_lclip_cnt>0 || $cns_ldisc_cnt>0 || $cns_rdisc_cnt>0){ #left part of a read was clipped, it's the right side of a L1
			$h_bk{left_part_clipped}{$chr}{$pos} = $lclip_cnt;
#		}
	}
}

#L1_pos  cns_lclip_cnt   cns_rclip_cnt   cnt_lpolyA_reads_from_sam_zx    cnt_rpolyA_reads_from_sam_zx    lpolyA_in_cns   rpolyA_in_cns   raw_ldisc_cnt   raw_rdisc_cnt   cns_ldisc_cnt   cns_rdisc_cnt   ref_repeat      cnv_pos l_distance      r_distance      cnv_len cnv_type
open CP, $candidate_pseudogene or die $!;
while(<CP>){
	chomp;
	if(/pos/){
		print "$_\n";
		next;
	}
	my @a = split;
	my $pos_info = $a[1];
	my ($chr, $pos1, $pos2) = $pos_info=~/(\w+):(\d+)-(\d+)/;
	my $flag1 = 0; #indicator of whether there are nearby clip site with the opposite clip direction
    my $flag2 = 0; #indicator of whether there are nearby clip site with the same clip direction
	my @tmp = (1, 2);
	my $pos;
	for my $index(@tmp){
		if(int($index) == 1){
			$pos = $pos1;
		}elsif(int($index) == 2){
			$pos = $pos2;
		}
###############not using########################################
		if(exists $h_bk{left_part_clipped}{$chr}{$pos}){
			for my $i($pos-50 .. $pos+50){
				if(exists $h_bk{right_part_clipped}{$chr}{$i} && $h_bk{right_part_clipped}{$chr}{$i}>=2){
					$flag1 = 1;
				}
			}
		}
		if(exists $h_bk{right_part_clipped}{$chr}{$pos}){
			for my $i($pos-50 .. $pos+50){
				if(exists $h_bk{left_part_clipped}{$chr}{$i} && $h_bk{left_part_clipped}{$chr}{$i}>=2){
					$flag1 = 1;
				}
			}
		}
#################################################################

		if(exists $h_bk{left_part_clipped}{$chr}{$pos}){
            for my $i($pos-10 .. $pos+10){
                if(exists $h_bk{left_part_clipped}{$chr}{$i} && $h_bk{left_part_clipped}{$chr}{$i}>=2){
					next if $i == $pos;
                    $flag2 = 1; #indicate that there are clipped reads with the same clip direction with 10bps of the target clip site
                }
            }
        }
        if(exists $h_bk{right_part_clipped}{$chr}{$pos}){
            for my $i($pos-10 .. $pos+10){
                if(exists $h_bk{right_part_clipped}{$chr}{$i} && $h_bk{right_part_clipped}{$chr}{$i}>=2){
					next if $i == $pos;
                    $flag2 = 1;
                }
            }
        }	
	}	
	
#	next if ($flag1==1); #next if there is nearby clip site with the opposite clip direction near the target cancidate clip site
	if($flag2 == 0){ #both 5' or 3' endpoints don't have nearby clip sites with the same clip direction
        print "$_\n";
    }else{
#		print "$_\n";
	}
		
}
close CP;	
