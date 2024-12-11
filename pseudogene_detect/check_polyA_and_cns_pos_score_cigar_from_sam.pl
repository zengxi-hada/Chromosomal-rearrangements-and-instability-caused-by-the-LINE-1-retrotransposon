#!/usr/bin/perl

use strict;
use warnings;

my $sam = shift; #L1_CL_GFPOrf_pos_A2.bam.clipped.sam
my $clip_sites = shift; #candidate_list_from_clip.txt_tmp

#######################################################################################################################################
## this program is to add the mapping position of supporting reads on the LINE1 cns											         ##
#######################################################################################################################################

die "perl $0 <*.bam.clipped.sam> <candidate_list_from_clip.txt_tmp>\n" unless ($sam && $clip_sites);

my %h_cns_pos; #hash containing polyA reads in L1 consensus sequence
my %h_cns_strand;
my %h_cns_cigar;
my %h_cns_polyA;
my %h_cns_algn_score;
my %h_cns_mismatch;
open SAM1, $sam or die $!;
while(<SAM1>){
	next if /^@/;
	chomp;
	my @a = split;
	my $read_id = $a[0];
	my $flag = $a[1];
	my $ref = $a[2];
	my $cns_pos = $a[3];
	my $align_score = $a[4];
	my $cigar = $a[5];
	my $seq = $a[9];
	my $cnt_mismatch = $a[11];
	$cnt_mismatch =~ s/NM:i://;
	next if $cigar =~ /H/;
	$h_cns_cigar{$read_id} = $cigar;
	$h_cns_algn_score{$read_id} = $align_score;
	$h_cns_mismatch{$read_id} = $cnt_mismatch;
	my $strand;
#	if($flag==0 || $flag==2048){
#		$strand = "positive";
#	}elsif($flag==16 || $flag==2064){
#		$strand = "negative";
#	}
	if($flag & 0x10){
        $strand = "negative";
    }else{
        $strand = "positive";
    }
	my $match_len = 0;
	while($cigar =~/(\d+)[DM]/g){
		$match_len += $1;
	}
#	$h_cns_pos{$read_id} = "$cns_pos\t$matched_bases";
	#polyA reads will not follow the rule below. because polyA has the same seq at both end, it will map to any pos to the ref polyA sequence
	if($read_id =~ /~L/){ #left clipped for original intact reads, need to collect the downstream pos as the joint pos
		if($strand eq "positive"){
			$h_cns_pos{$read_id} = $cns_pos + $match_len;
			$h_cns_strand{$read_id} = "+";
		}else{
			$h_cns_pos{$read_id} = $cns_pos;
			$h_cns_strand{$read_id} = "-";
		}
	}elsif($read_id =~ /~R/){ #right clipped for original intact reads, need to collect the upstream pos as the joint pos
		if($strand eq "positive"){
            $h_cns_pos{$read_id} = $cns_pos;
			$h_cns_strand{$read_id} = "+";
        }else{
            $h_cns_pos{$read_id} = $cns_pos + $match_len;
			$h_cns_strand{$read_id} = "-";
        }
	}
	
	if($flag & 0x10){ #negative strand
		$seq = reverse($seq);
		$seq =~ tr/ATCGatcg/TAGCtagc/;
	}
	my $len = length($seq); #
	my $cnt_A = $seq =~ tr/Aa/Aa/;
	my $cnt_T = $seq =~ tr/Tt/Tt/;
	my $A_ratio = $cnt_A / $len;
	my $T_ratio = $cnt_T / $len;
	if($A_ratio>0.8){
		$h_cns_polyA{polyA}{$read_id} = 1;
	}elsif($T_ratio>0.8){
		$h_cns_polyA{polyT}{$read_id} = 1;
	}
	
#	if($ref =~ /_polyA/){
#		$h_cns_polyA{$read_id} = 1; # the clip part of the original intact read is polyA
#	}	
}
close SAM1;

print "chr\tpos\tlclip_cnt\trclip_cnt\tcnt_mate_in_cns\tcns_lclip_cnt\tcns_rclip_cnt\tcns_lpos,cns_rpos\n";
open CS, $clip_sites or die $!;
while(<CS>){
	chomp;
	my @a = split;
    my $chr = $a[0];
    my $pos = $a[1];
    my $lclip_cnt = $a[2]; #raw count
    my $rclip_cnt = $a[3]; #raw count
#	next if ($lclip_cnt+$rclip_cnt) < 2;
    my $cnt_mate_in_cns = $a[4];
    my $cns_lclip_cnt = $a[5];# reads whose left clipped part could map to LINE1-cns
    my $cns_rclip_cnt = $a[6];# reads whose right clipped part could map to LINE1-cns
	my $cns_reads_str = $a[8]; # LINE-1 consensus clipped reads
	my $general_reads_str = $a[7]; # general clipped reads
	my @b = split /;/, $cns_reads_str; # split left and right clipped reads
	my (@left_cns_algn_score, @right_cns_algn_score, @left_cns_cigar, @right_cns_cigar, @left_cns_polyA, @right_cns_polyA, @left_cns_pos, @right_cns_pos, @left_cns_strand, @right_cns_strand);
	my (@left_cns_mismatch, @right_cns_mismatch);
	if($cns_reads_str =~ /\w+;/){ # reads whose left clipped part could map to LINE1-cns, 3' end of the L1 insertion
		my @left_reads = split /,/, $b[0]; # split reads within left or within right side
		for my $i(@left_reads){ #reads whose left clipped part could map to LINE1-cns, right side of the LINE1 insertion
			my $cns_pos = $h_cns_pos{$i};
			my $strand = $h_cns_strand{$i};
			my $cigar = $h_cns_cigar{$i};
			my $align_score = $h_cns_algn_score{$i};
			my $cnt_mismatch = $h_cns_mismatch{$i};
			my $polyA_or_not = (exists $h_cns_polyA{polyA}{$i} || exists $h_cns_polyA{polyT}{$i}) ? "polyA" : "not_polyA";
#			my ($cns_pos, $matched_bases) = split /\t/, $cns_info;
			push @left_cns_pos, $cns_pos;
			push @left_cns_strand, $strand;
			push @left_cns_algn_score, $align_score;
			push @left_cns_cigar, $cigar;
			push @left_cns_polyA, $polyA_or_not;
			push @left_cns_mismatch, $cnt_mismatch;
#			push @left_cns_match_num, $matched_bases;
			
		}
	}elsif($cns_reads_str =~ /;\w+/){ #reads whose right clipped part could map to LINE1-cns, 5' end of the L1 insertion
		my @right_reads = split /,/, $b[1]; # split reads within left or within right side
		for my $i(@right_reads){ #reads whose right clipped part could map to LINE1-cns, left side of the LINE1 insertion
			my $cns_pos = $h_cns_pos{$i};
			my $strand = $h_cns_strand{$i};
			my $cigar = $h_cns_cigar{$i};
			my $align_score = $h_cns_algn_score{$i};
			my $cnt_mismatch = $h_cns_mismatch{$i};
			my $polyA_or_not = (exists $h_cns_polyA{polyA}{$i} || exists $h_cns_polyA{polyT}{$i}) ? "polyA" : "not_polyA";
#			my ($cns_pos, $matched_bases) = split /\t/, $cns_info;
			push @right_cns_pos, $cns_pos;
			push @right_cns_strand, $strand;
			push @right_cns_algn_score, $align_score;
			push @right_cns_cigar, $cigar;
			push @right_cns_polyA, $polyA_or_not;
			push @right_cns_mismatch, $cnt_mismatch;
#			push @right_cns_match_num, $matched_bases;
		}
	}
#	@left_cns_pos = (0) if(not @left_cns_pos); #no supporting reads of this site can map to LINE1
#	@right_cns_pos = (0) if(not @right_cns_pos); #no supporting reads of this site can map to LINE1
#	my $mean_left_cns_pos = &mean(@left_cns_pos);
#	my $mean_right_cns_pos = &mean(@right_cns_pos);
	my $most_frq_left_cns_pos = &get_the_most_frq_item(@left_cns_pos);
	my $most_frq_right_cns_pos = &get_the_most_frq_item(@right_cns_pos);
#	my $mean_left_match_num = &mean(@left_cns_match_num);
#   my $mean_right_match_num = &mean(@right_cns_match_num);
	my $left_cns_pos_str = join(",", @left_cns_pos);
	my $right_cns_pos_str = join(",", @right_cns_pos);
	my $left_cns_strand_str = join(",", @left_cns_strand);
    my $right_cns_strand_str = join(",", @right_cns_strand);
	my $left_cns_algn_score_str = join(",", @left_cns_algn_score);
    my $right_cns_algn_score_str = join(",", @right_cns_algn_score);
	my $left_cns_cigar_str = join(",", @left_cns_cigar);
    my $right_cns_cigar_str = join(",", @right_cns_cigar);
	my $left_cns_polyA_str = join(",", @left_cns_polyA);
    my $right_cns_polyA_str = join(",", @right_cns_polyA);
	my $left_cns_mismatch_str = join(",", @left_cns_mismatch);
	my $right_cns_mismatch_str = join(",", @right_cns_mismatch);
#	print "$chr\t$pos\t$lclip_cnt\t$rclip_cnt\t$cnt_mate_in_cns\t$cns_lclip_cnt\t$cns_rclip_cnt\t$refined_polyA_or_T_reads\n";
	print "$chr\t$pos\t$lclip_cnt\t$rclip_cnt\t$cnt_mate_in_cns\t$cns_lclip_cnt\t$cns_rclip_cnt\t$most_frq_left_cns_pos,$most_frq_right_cns_pos\t$left_cns_pos_str;$right_cns_pos_str\t$left_cns_strand_str;$right_cns_strand_str\t$left_cns_algn_score_str;$right_cns_algn_score_str\t$left_cns_cigar_str;$right_cns_cigar_str\t$left_cns_mismatch_str;$right_cns_mismatch_str\t$left_cns_polyA_str;$right_cns_polyA_str\t$general_reads_str\n";
}
close CS;

################# subroutine #####################
sub mean
{
    my @a = @_;
    my $ave; my $sum;
	my $count = 0;
	if(not @a){
		return "nan";
	}
    for my $i(@a){
   		$sum += $i;
		$count++;
	}
	$ave = int($sum / $count);
    return $ave;
}

sub get_the_most_frq_item
{
    my @a = @_;
    my $count = 0;
    my %h;
	if(not @a){
        return "nan";
    }
    for my $item(@a){
#        $count++;
        $h{$item}++;
    }
	my %frq;
	for my $item(keys %h){
		my $count = $h{$item};
		$frq{$count} = $item;
	}
    my $most_frq_item;
    for my $count(sort {$b<=>$a} keys %frq){
        $most_frq_item = $frq{$count};
        last;
    }
    return $most_frq_item;
}
