#!/usr/bin/perl

use strict;
use warnings;

my $clip_disc_info = shift; #clip_disc_polyA_Info.txt.exclude_L1_info

die "perl $0 <clip_disc_polyA_info.txt>\n" unless ($clip_disc_info);

my %h_record;
open CD, $clip_disc_info or die $!;
while(<CD>){
	if(/nlclip/ || /nlcip/){
		print "chr\tpos\tnlclip\tnrclip\tnldisc\tnrdisc\tpolyA_info\tmost_frq_supple_algn\tmost_frq_mate_of_clip_reads\tmost_frq_mate_of_disc_reads\tsupple_algn_str\tmate_of_clip_reads_str\tmate_of_disc_reads_str\n";
		next;
	}
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $pos = $a[1];
	my $nlclip = $a[2];
	my $nrclip = $a[3];
	my $nldisc = $a[4];
	my $nrdisc = $a[5];
	my $polyA_info = $a[6];
	my $supple_algn = $a[7];
	my $clip_reads = $a[8];
	my $disc_reads = $a[9];
	if($disc_reads eq ";"){
		$disc_reads = "nan;nan";
	}elsif($disc_reads =~ /^;/){
		$disc_reads = "nan".$disc_reads;
	}elsif($disc_reads =~ /;$/){
		$disc_reads = $disc_reads."nan";
	}
	my @supple_algn_arr_tmp = split /,/, $supple_algn;
	my @clip_reads_arr = split /,/, $clip_reads;
    my @disc_reads_arr1 = split /;/, $disc_reads;
	my @disc_reads_arr2 = ($nlclip==0) ? (split /,/, $disc_reads_arr1[0]) : (split /,/, $disc_reads_arr1[1]);
	my (@supple_algn_arr, @mate_of_clip_reads, @mate_of_disc_reads);
	for my $i(@supple_algn_arr_tmp){
		next if $i =~ /nan/;
		push @supple_algn_arr, $i;
	}

	for my $i(@clip_reads_arr){
#        next if $i=~/nan$/;
		my @tmp = split /~/, $i;
		my $mate_chr = $tmp[3];
		my $mate_pos = $tmp[4];
		push @mate_of_clip_reads, "$mate_chr:$mate_pos";
	}
	
	for my $i(@disc_reads_arr2){
		next if $i eq "nan";
		my @tmp = split /~/, $i;
        my $mate_chr = $tmp[2];
        my $mate_pos = $tmp[3];
        push @mate_of_disc_reads, "$mate_chr:$mate_pos";
	}

	@mate_of_disc_reads=("no_disc") if @mate_of_disc_reads==0;
	@supple_algn_arr=("no_supple") if @supple_algn_arr==0;	

	my $supple_algn_str = join (",", @supple_algn_arr);
	my $mate_of_clip_reads_str = join (",", @mate_of_clip_reads);
	my $mate_of_disc_reads_str = join (",", @mate_of_disc_reads);

	my ($most_frq_supple_algn1, $most_frq_supple_algn2) = &get_two_most_frq_item(@supple_algn_arr);
	my ($most_frq_mate_of_clip_reads1, $most_frq_mate_of_clip_reads2) = &get_two_most_frq_item(@mate_of_clip_reads);
	my ($most_frq_mate_of_disc_reads1, $most_frq_mate_of_disc_reads2) = &get_two_most_frq_item(@mate_of_disc_reads);	
	print "$chr\t$pos\t$nlclip\t$nrclip\t$nldisc\t$nrdisc\t$polyA_info\t$most_frq_supple_algn1,$most_frq_supple_algn2\t$most_frq_mate_of_clip_reads1,$most_frq_mate_of_clip_reads2\t$most_frq_mate_of_disc_reads1,$most_frq_mate_of_disc_reads2\t$supple_algn_str\t$mate_of_clip_reads_str\t$mate_of_disc_reads_str\n";
	
#	$h_record{$chr}{$pos} = "$most_frq_supple_algn1\t$most_frq_supple_algn2\t$most_frq_mate_of_clip_reads1\t$most_frq_mate_of_clip_reads2\t$most_frq_mate_of_disc_reads1\t$most_frq_mate_of_disc_reads2";
#@	$h_record{$chr}{$pos} = "$most_frq_supple_algn1\t$mate_of_disc_reads_str";
	
}


################# subroutine ##############
sub get_two_most_frq_item
{
    my @a = @_;
    my $count = 0;
    my %h;
    if(not @a){
        return "nan";
    }
    for my $item(@a){
        $h{$item}++;
    }
    my %frq;
    for my $item(keys %h){
        my $count = $h{$item};
        $frq{$count} = $item;
    }
    my $most_frq_item;
	my $count2 = 0;
	my %h2;
    for my $count(sort {$b<=>$a} keys %frq){
		$count2++;
        $most_frq_item = $frq{$count};
		$h2{$count2} = $most_frq_item;
        last if $count2 == 2;;
    }
	if($count2 == 1){
	    return ($h2{1}, $h2{1});
	}else{
		return ($h2{1}, $h2{2});
	}
}

sub get_two_most_frq_disc_item
{
    my @a = @_;
    my $count = 0;
    my %h;
    if(not @a){
        return "nan";
    }
    for my $item(@a){
        $h{$item}++;
    }
    my %frq;
	my %count_new;
    for my $item(keys %h){
		next if exists $count_new{$item};
		for my $item_tmp(keys %h){
			next if $item eq $item_tmp;
			next if exists $count_new{$item_tmp};
			my ($chr, $pos) = split /:/, $item;
			my ($chr_tmp, $pos_tmp) = split /:/, $item_tmp;
			if(($chr eq $chr_tmp) && abs($pos-$pos_tmp)<10000){
				$count_new{$item} += 2; #the number of occurence for a pos
			}
		}
    }

	for my $item(keys %count_new){
		my $count = $count_new{$item};
        $frq{$count} = $item;
	}	

    my $most_frq_item;
    my $count2 = 0;
    my %h2;
    for my $count(sort {$b<=>$a} keys %frq){
        $count2++;
        $most_frq_item = $frq{$count};
        $h2{$count2} = $most_frq_item;
        last if $count2 == 2;;
    }
    if($count2 == 1){
        return ($h2{1}, $h2{1});
    }else{
        return ($h2{1}, $h2{2});
    }
}
