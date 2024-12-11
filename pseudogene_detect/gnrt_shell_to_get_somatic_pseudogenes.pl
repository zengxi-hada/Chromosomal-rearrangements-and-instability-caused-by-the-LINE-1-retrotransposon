#!/usr/bin/perl

use strict;
use warnings;

my $file_list = shift; #vcf list for sideRetro

die "perl $0 <file.list>\n" unless $file_list;

open FL, $file_list or die $!;
while(<FL>){
	chomp;
	my $vcf = $_;
	my @a = split /\//, $vcf;
	my $sample = $a[-2];
	$sample =~ s/-sideretro//;
	open OUT, ">$sample/L1/get_somatic_pseudogene.sh" or die $!;
    print OUT "#!/bin/bash\n";
    print OUT "#SBATCH -c 2\n";
    print OUT "#SBATCH -t 0-5:10\n";
    print OUT "#SBATCH --mem=60G\n";
    print OUT "#SBATCH -p bch-compute\n";
    print OUT "#SBATCH -o $sample\_get_somatic_%j.out\n\n";

#@	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/flt_sideRetro_set_v2.pl $vcf /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.final > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.final.intesect_with_sideretro\n";
#@	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/flt_sideRetro_set_v2.pl $vcf /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.intesect_with_sideretro\n";
#@	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/flt_sideRetro_set_v2.pl $vcf /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.intesect_with_sideretro\n";
	
#@	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/get_somatic_pseudogene.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/all_DMSO.list /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.final.intesect_with_sideretro > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_somatic_v1.txt\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/get_somatic_pseudogene.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/all_clonal_DMSO.list /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.intesect_with_sideretro > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_somatic_v2.txt\n";
#@	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/get_somatic_pseudogene.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/all_DMSO.list /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_rm_rep.txt.intesect_with_sideretro > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_somatic_v5.txt\n";
#@	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/get_somatic_pseudogene.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/all_DMSO.list /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes.txt > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/$sample/L1/pseudogenes_somatic_v6.txt\n";

}
close FL;
close OUT;
