#!/usr/bin/perl

use strict;
use warnings;

my $full_clip_file = shift; #candidate_list_from_clip.txt_tmp.add_cns_pos
my $full_disc_file = shift; #candidate_list_from_disc.txt 
my $query_file = shift; #ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio, where some clip pos has been merged into one single breakpoint (peak)

die "perl $0 <candidate_list_from_clip.txt_tmp.add_cns_pos> <candidate_list_from_disc.txt> <ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio>\n" unless ($full_clip_file && $full_disc_file && $query_file);

my %h_clip;
open IN1, $full_clip_file or die $!;
while(<IN1>){
	next if /lclip_cnt/;
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $LINE1_inside_pos = $a[7];
	my $LINE1_inside_pos_str = $a[8];
	my $clip_part_strand = $a[9];
	my $clip_part_LINE1_algn_score = $a[10];
	my $clip_part_LINE1_algn_cigar = $a[11];
	my $clip_part_LINE1_algn_mismatch = $a[12];
	my $polyA_info = ($a[13] eq ";") ? "not_polyA" : $a[13];
	my $reads_str = $a[14]; #general reads str
	$reads_str =~ s/;/,/;
	my @reads_arr = split /,/, $reads_str;
	my @supple_algn_arr = ();
	for my $read(@reads_arr){
        my $origin_read = $read;
		my @b = split /~/, $read;
        my $supple_algn_chr = $b[-2];
		my $supple_algn_pos = $b[-1];
        $read = "$b[0]~$b[1]";
		my $supple_algn = ($origin_read=~/nan$/) ? "nan" : "$supple_algn_chr:$supple_algn_pos";
		push @supple_algn_arr, $supple_algn;
	}
	my $supple_algn_str = join (",", @supple_algn_arr); 

	$h_clip{"$chr\t$pos"} = "$LINE1_inside_pos\t$LINE1_inside_pos_str\t$clip_part_strand\t$clip_part_LINE1_algn_score\t$clip_part_LINE1_algn_cigar\t$polyA_info\t$supple_algn_str\t$reads_str";
#	$h_clip{"$chr\t$pos"} = "$clip_part_LINE1_algn_score\t$clip_part_LINE1_algn_cigar\t$polyA_info\t$reads_str";
}
close IN1;

my %h_disc; #it is actually the intersection of disc and clip
open FDF, $full_disc_file or die $!;
while(<FDF>){
	chomp;
	my @a = split;
	my $chr = $a[0];
    my $pos = $a[1];
	my $disc_reads = $a[-1];
	$h_disc{"$chr\t$pos"} = $disc_reads;
}
close FDF;

print "type\tchr\tpos\tnlclip\tnrclip\tn_mate_in_L1_cns\tL1_cns_nlclip\tL1_cns_nrclip\tLINE1_ratio_of_clip_reads\tpolyA_info\tlpolyA_in_cns\trpolyA_in_cns\tn_ldisc\tn_rdisc\tn_L1_cns_ldisc\tn_L1_cns_rdisc\tL1_inside_pos\tL1_inside_pos_str\tclip_part_L1_strand\tclip_part_LINE1_algn_score\tclip_part_LINE1_algn_cigar\tpolyA_info\tsupple_algn_of_clip_reads\tclip_reads\tdisc_reads\n";
open QF, $query_file or die $!;
while(<QF>){
	next if /lclip_cnt/;
	chomp;
	my @a = split;
	my $chr = $a[1];
    my $pos = $a[2];
	my $clip_info = $h_clip{"$chr\t$pos"};
	my $disc_info = (exists $h_disc{"$chr\t$pos"}) ? $h_disc{"$chr\t$pos"} : "nan";
	print "$_\t$clip_info\t$disc_info\n";
}
close QF;
