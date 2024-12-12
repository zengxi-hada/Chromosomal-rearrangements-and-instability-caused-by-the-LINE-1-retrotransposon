#!/usr/bin/perl

use strict;
use warnings;

my $in = shift; #ORFeus_insertion.2-types.txt.add_clip_ORFeus_ratio.flt.refine
my $vcf = shift; #*LINE1.vcf

die "perl $0 <candidate.pre> <*LINE1.vcf>\n" unless ($in && $vcf);

my %candidate;
my %tmp;
open IN, $in or die $!;
while(<IN>){
	next if /filtered/;
    chomp;
    my @a = split;
    my $pos_str = $a[0];
    my ($chr, $lpos, $rpos) = $pos_str=~ /(\w+):(\d+)-(\d+)/;
	$candidate{$chr}{$lpos} = $rpos;
	$tmp{$chr}{$lpos} = $rpos;
    $tmp{$chr}{$rpos} = $lpos;
}
close IN;

open VCF, $vcf or die $!;
while(<VCF>){
	chomp;
	next if /^#/;
    chomp;
	my @a = split;
	my $chr = $a[0];
	my $lpos = $a[1];
	my $info = $a[7];
	my @tmp = split /;/, $info;
	my $rpos = $tmp[2];
	$rpos =~ s/END=//;
	if((not exists $tmp{$chr}{$lpos}) && (not exists $tmp{$chr}{$rpos})){
		$candidate{$chr}{$lpos} = "$rpos#";
	}
}
close VCF;

for my $chr(keys %candidate){
	for my $lpos(keys %{$candidate{$chr}}){
		my $rpos = $candidate{$chr}{$lpos};
		$rpos =~ s/#//;
		print "$chr\t$lpos\t$rpos\n";
	}
}
