#!/usr/bin/perl

use strict;
use warnings;

###########################################################################################################################################
# This program will filter out some unqualified candidate breakpoints for pseudogene insertion, then extract the psedugene ins.          ##
# The clipped part of clip supporting reads and mate of discordant supporting reads should map to a close location on the human genome.  ##
###########################################################################################################################################

my $clip_disc_bk = shift; #clip_disc_polyA_Info.txt.exclude_L1_info.ts.sort

die "perl $0 <clip_disc_polyA_Info.txt.exclude_L1_info.ts>\n" unless $clip_disc_bk;

#get the qualified breakpoints for pseudogene
my %h_bk; # store the breakpoints
open CDB, $clip_disc_bk or die $!;
while(<CDB>){
	next if /nrclip/;
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $nlclip = $a[2];
	my $nrclip = $a[3];
	my $nldisc = $a[4];
	my $nrdisc = $a[5];
	my $polyA_info = $a[6];
	$polyA_info =~ s/;//g;
	my $two_most_frq_supple = $a[7];
	my $all_supple = $a[10];
	my $clip_mates = $a[11];
	my $disc_mates = $a[12];
	my @two_most_frq_supple_arr = split /,/, $two_most_frq_supple;
	my ($most_frq_supple_algn_chr, $most_frq_supple_algn_pos) = split /:/, $two_most_frq_supple_arr[0];
	my @disc_mates_arr = split /,/, $disc_mates;
	my $count_hits = 0;
	for my $i(@disc_mates_arr){
		my ($disc_mate_chr, $disc_mate_pos) = split /:/, $i;
		#the supplementary alignment pos and mate reads align pos should be closei (within 20000 bps)
		if(($most_frq_supple_algn_chr eq $disc_mate_chr) && abs($most_frq_supple_algn_pos-$disc_mate_pos)<20000){
			$count_hits++;
		} 
	}
	if($polyA_info !~ /not_polyA/){#when it is a polyA breakpoint, then no need to filter, just use it
		if($nrclip > 0){ #bk supported by right clipped reads
            $h_bk{right_clip}{"$chr\t$pos"} = "$nrclip\t$nldisc\t$polyA_info\t$two_most_frq_supple\t$all_supple\t$clip_mates\t$disc_mates";
        }elsif($nlclip > 0){ #bk supported by left clipped reads
            $h_bk{left_clip}{"$chr\t$pos"} = "$nlclip\t$nrdisc\t$polyA_info\t$two_most_frq_supple\t$all_supple\t$clip_mates\t$disc_mates";
        }
	}elsif($polyA_info=~/not_polyA/ && $count_hits >= 1){ #when not_polyA, the supplementary alignment pos and mate reads align pos are close, a signal for pseudogene insertion breakpoint
#		print "$_\n";
		if($nrclip > 0){ #bk supported by right clipped reads
	        $h_bk{right_clip}{"$chr\t$pos"} = "$nrclip\t$nldisc\t$polyA_info\t$two_most_frq_supple\t$all_supple\t$clip_mates\t$disc_mates";
	    }elsif($nlclip > 0){ #bk supported by left clipped reads
	        $h_bk{left_clip}{"$chr\t$pos"} = "$nlclip\t$nrdisc\t$polyA_info\t$two_most_frq_supple\t$all_supple\t$clip_mates\t$disc_mates";
	    }
	}
}
close CDB;

#pair the breakpoints
print "sides\tchr:pos\tclip_cnt\tdisc_cnt\tpolyA_info\ttwo_most_frq_supple\tall_supple_algn\tclip_mates\tdisc_mates\n";
my %left_only; # store the pseudogene insertions with only left side signal
#print for pseudogene insertions with both left and right clip supporting reads
for my $lpos_str(keys %{$h_bk{left_clip}}){
    my $flag = 0;
    my @l_arr = split /\t/, $lpos_str;
    my $lchr = $l_arr[0];
    my $lpos = $l_arr[1];
    my @l_cnt_arr = split /\t/, $h_bk{left_clip}{$lpos_str};
    my $lclip_cnt = $l_cnt_arr[0];
    my $lpolyA_info = $l_cnt_arr[2];
    my $ldisc_cnt = $l_cnt_arr[1];
	my $l_two_most_frq_supple = $l_cnt_arr[3];
	my $l_supple = $l_cnt_arr[4];
	my $l_clip_mates = $l_cnt_arr[5];
	my $l_disc_mates = $l_cnt_arr[6];
    for my $rpos_str(keys %{$h_bk{right_clip}}){
        my @r_arr = split /\t/, $rpos_str;
        my $rchr = $r_arr[0];
        my $rpos = $r_arr[1];
        if(($lchr eq $rchr) && abs($lpos-$rpos)<=200){ #the left and right breakpoint should be <= 200 in a pseudogene insertion
            my @r_cnt_arr = split /\t/, $h_bk{right_clip}{$rpos_str};
            my $rclip_cnt = $r_cnt_arr[0];
            my $rpolyA_info = $r_cnt_arr[2];
            my $rdisc_cnt = $r_cnt_arr[1];
			my $r_two_most_frq_supple = $r_cnt_arr[3];
			my $r_supple = $r_cnt_arr[4];
		    my $r_clip_mates = $r_cnt_arr[5];
		    my $r_disc_mates = $r_cnt_arr[6];
            print "left-right\t$lchr:$lpos-$rpos\t$lclip_cnt,$rclip_cnt\t$ldisc_cnt,$rdisc_cnt\t$lpolyA_info;$rpolyA_info\t$l_two_most_frq_supple;$r_two_most_frq_supple\t$l_supple;$r_supple\t$l_clip_mates;$r_clip_mates\t$l_disc_mates;$r_disc_mates\n";
            $flag = 1;
        }
    }
    if($flag==0 && $lclip_cnt>=1){ #when not right clip breakpoint can pair with the current left clip bk
    	$left_only{$lpos_str} = $h_bk{left_clip}{$lpos_str};
    }                  
}

##for the insertions only have right or left side
#for breakpoint with both ORFeus clip and disc reads
my %right_only;
for my $rpos_str(keys %{$h_bk{right_clip}}){
	my $flag = 0;
	my @r_arr = split /\t/, $rpos_str;
	my $rchr = $r_arr[0];
	my $rpos = $r_arr[1];
	my @r_cnt_arr = split /\t/, $h_bk{right_clip}{$rpos_str};
    my $rclip_cnt = $r_cnt_arr[0];
    my $rpolyA_info = $r_cnt_arr[2];
    my $rdisc_cnt = $r_cnt_arr[1];
    for my $lpos_str(keys %{$h_bk{left_clip}}){
		my @l_arr = split /\t/, $lpos_str;
        my $lchr = $l_arr[0];
        my $lpos = $l_arr[1];
		if(($lchr eq $rchr) && abs($lpos-$rpos)<=200){
			$flag = 1;
		}
	}
	if($flag==0 && $rclip_cnt>=1){ #when no left clip pos could map to the current right clip pos, then store it as right only bk
		$right_only{$rpos_str} = $h_bk{right_clip}{$rpos_str};
	}
}

#print pseudogene insertions with only left breakpoint
for my $pos_str(keys %left_only){
	my @a = split /\t/, $pos_str;
	my $chr = $a[0];
	my $pos = $a[1];
	my $info = $left_only{$pos_str};
	print "only_left\t$chr:$pos-$pos\t$info\n";
}

#print pseudogene insertions with only right breakpoint
for my $pos_str(keys %right_only){
    my @a = split /\t/, $pos_str;
    my $chr = $a[0];
    my $pos = $a[1];
	my $info = $right_only{$pos_str};
    print "only_right\t$chr:$pos-$pos\t$info\n";
}
