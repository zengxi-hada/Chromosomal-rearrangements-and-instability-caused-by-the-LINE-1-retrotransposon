#!/usr/bin/perl

use strict;
use warnings;

my $sam1 = shift; #polyA sam, *.sam_partial_polyA.sam
my $sam2 = shift; #fully mapped sam of reads whoses clipped part could fully mapped to cns, *.bam.clipped.sam.sam
my $clip_sites = shift; #candidate_list_from_clip.txt_tmp

#####################################################################################################################################################################
## this program is to refine the polyA info, which is to further add polyA info directly from sam file of realignment of cns and LIHS(repeat) region in clip step  ##
#####################################################################################################################################################################

die "perl $0 <.sam_partial_polyA.sam> <clipped_part_fully_map.sam> <candidate_list_from_clip.txt_tmp>\n" unless ($sam1 && $sam2 && $clip_sites);

my %h; #hash containing polyA reads
open SAM1, $sam1 or die $!;
while(<SAM1>){
	next if /^@/;
	chomp;
	my @a = split;
	my $read_id = $a[0];
	$h{$read_id} = 1;
}
close SAM1;

open SAM2, $sam2 or die $!;
while(<SAM2>){
    next if /^@/;
    chomp;
    my @a = split;
    my $read_id = $a[0];
	my $flag = $a[1];
	my $seq = $a[9];
#	my $strand = $a[1]; #0, positive strand; 16, negative strand
#	if($strand == 16){	#negative strand
	if($flag & 0x10){ #negative strand
		$seq = reverse($seq);
		$seq =~ tr/ATCGatcg/TAGCtagc/;
	}
	my $len = length($seq); #
	my $cnt_A = $seq =~ tr/Aa/Aa/;
	my $cnt_T = $seq =~ tr/Tt/Tt/;
	my $A_ratio = $cnt_A / $len;
	my $T_ratio = $cnt_T / $len;
#	print "##$read_id\t$seq\t$cnt_A\t$cnt_T\t$A_ratio\t$T_ratio\n";
	if($A_ratio>0.8){
		$h{polyA}{$read_id} = 1;
#		print "##$read_id\t$seq\t$cnt_A\t$cnt_T\t$A_ratio\t$T_ratio\n";
	}elsif($T_ratio>0.8){
		$h{polyT}{$read_id} = 1;
	}
}
close SAM1;

print "chr\tpos\tlclip_cnt\trclip_cnt\tcnt_mate_in_cns\tcns_lclip_cnt\tcns_rclip_cnt\tcnt_polyA_reads_from_sam_zx\n";
open CS, $clip_sites or die $!;
while(<CS>){
	chomp;
	my @a = split;
    my $chr = $a[0];
    my $pos = $a[1];
    my $lclip_cnt = $a[2]; #raw count
    my $rclip_cnt = $a[3]; #raw count
#	next if ($lclip_cnt+$rclip_cnt) < 2; #only consider the breakpoint with at least 2 supporting reads
    my $cnt_mate_in_cns = $a[4];
    my $cns_lclip_cnt = $a[5];# reads whose left clipped part could map to LINE1-cns
    my $cns_rclip_cnt = $a[6];# reads whose right clipped part could map to LINE1-cns
	my $reads_str = $a[7]; # clipped reads whose clipped part could map to cns(LINE1, LNE1-2 ...) or LIHS (repeat region)
#	my $refined_polyA_or_T_reads = 0;
	my $refined_polyA_reads = 0;
	my $refined_polyT_reads = 0;
	my @b = split /;/, $reads_str; # split left and right clipped reads
	if($reads_str =~ /\w+;/){ # reads whose left clipped part could map to LINE1-cns, 3' end of the L1 insertion
		my @left_reads = split /,/, $b[0]; # split reads within left or within right side
		for my $i(@left_reads){ #reads whose left clipped part could map to LINE1-cns, right side of the LINE1 insertion
			if(exists $h{polyA}{$i}){ #only polyA could happen when it is a real LINE1 insertion
				$refined_polyA_reads++;
			}elsif(exists $h{polyT}{$i}){
				$refined_polyT_reads++;
			}
		}
	}elsif($reads_str =~ /;\w+/){ #reads whose right clipped part could map to LINE1-cns, 5' end of the L1 insertion
		my @right_reads = split /,/, $b[1]; # split reads within left or within right side
		for my $i(@right_reads){ #reads whose right clipped part could map to LINE1-cns, left side of the LINE1 insertion
			if(exists $h{polyA}{$i}){
                $refined_polyA_reads++;
            }elsif(exists $h{polyT}{$i}){ #only polyT could happen when it is a real LINE1 insertion
                $refined_polyT_reads++;
            }
		}
	}
#	print "$chr\t$pos\t$lclip_cnt\t$rclip_cnt\t$cnt_mate_in_cns\t$cns_lclip_cnt\t$cns_rclip_cnt\t$refined_polyA_or_T_reads\n";
	print "$chr\t$pos\t$lclip_cnt\t$rclip_cnt\t$cnt_mate_in_cns\t$cns_lclip_cnt\t$cns_rclip_cnt\tpolyA,$refined_polyA_reads;polyT,$refined_polyT_reads\n";
}
close CS;



