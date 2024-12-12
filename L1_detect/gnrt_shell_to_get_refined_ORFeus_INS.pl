#!/usr/bin/perl

use strict;
use warnings;

#################################################################################################################################################################
## This program is to generate shell to obtain the refined ORFeus L1 insertion sites only using the output of clip and disc steps of xTea (1st and 2nd step)   ##
## This program does not check the long range ORFeus L1 insertion, which will be addressed when detecting L1-medidated SV                                      ##
#################################################################################################################################################################

my $bam_list = shift;

die "perl $0 <bam.list>\n" unless $bam_list;

#my $dir = "/lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/ORFeus_GFP_L1_ins/run_with_only_ORFeus_L1_as_cns/use_mRNA_rm_intron_and_puromycin";
#my $dir = "/lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/ORFeus_GFP_L1_ins/using_both_plasmid_and_mRNA";
#my $dir = "/lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/ORFeus_GFP_L1_ins/use_plasmid_rm_puromycin";
#my $dir = "/lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/ORFeus_GFP_L1_ins/use_mRNA_rm_intron_and_puromycin";
#my $dir = "/lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/ORFeus_Dox_L1_ins/backup/no_threshold_in_clip_and_disc_stp";
#my $dir = "/lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/ORFeus_Dox_L1_ins/backup/batch2";
#@my $all_non_ref_list = "DMSO_DoX.candidate_disc_filtered_cns.txt.list"; #xTea germline list
my $all_non_ref_list = "DMSO_DoX.candidate_disc_filtered_cns.txt.list.ts";
#@my $all_non_ref_list = "GFP.candidate_disc_filtered_cns.txt.list";
my $zx_germline_list = "ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.list";

open BL, $bam_list or die $!;
while(<BL>){
#	next if !/Dox/;
	chomp;
	my @a = split;
	my $sample = $a[0];
#	open OUT, ">$sample/L1/get_singleton_and_somatic_ORFeus_L1.sh" or die $!;
	open OUT, ">$sample/L1/get_singleton_and_somatic_ORFeus_L1.sh" or die $!;
#	open OUT, ">$sample/L1/label_for_somatic_L1.sh" or die $!;
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH -c 1\n";
	print OUT "#SBATCH -t 0-6:00\n";
	print OUT "#SBATCH --mem=25G\n";
	print OUT "#SBATCH -p bch-compute\n";
	print OUT "#SBATCH -o $sample\_refine_%j.out\n\n";
    
	print OUT "\n#get the somatic ORFeus insertions\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/check_polyA_from_sam.pl $sample/L1/tmp/$sample.bam.clipped.sam_partial_polyA.sam $sample/L1/tmp/$sample.bam.clipped.sam.sam $sample/L1/candidate_list_from_clip.txt_tmp > $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/add_disc_info_at_candidate_list_from_clip.txt_tmp.pl $sample/L1/tmp/raw_discordant_reads_tmp0_seplf $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam > $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/add_polyA_info_again.pl $sample/L1_add_polyA_at_clip_stp/candidate_list_from_clip.txt_tmp $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc > $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc.add_xTea_polyA\n";
	print OUT "perl -lane \'next if /merged_in_peak/; if(/lclip_cnt/){print \"type1/2\\t\$_\"; next}; if((\$F[5]+\$F[6])>0 && (\$F[12]+\$F[13])>0){print \"with_both_clip_disc\\t\$_\"}elsif(\$F[5]+\$F[6]>0){print \"with_only_clip\\t\$_\";}elsif(\$F[12]+\$F[13]>0){print \"with_only_disc\\t\$_\";}else{print \"other\\t\$_\"}\' $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc.add_xTea_polyA | sort > $sample/L1/ORFeus_insertion.sub-types.txt\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/cal_realign_cnt_2-types.pl $sample/L1/ORFeus_insertion.sub-types.txt $sample/L1/tmp/$sample.bam.clipped.sam > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio\n";
#@	print OUT "perl -lane \'if(/lclip_cnt/){print;next} my \@a=split /,/, \$F[8]; my \$ratio=\$a[2]; if((\$F[0] eq \"with_both_clip_disc\" && \$ratio>=0.5) || (\$F[0] eq \"with_only_disc\") || ((\$F[0] eq \"with_only_clip\") && (\$F[6]+\$F[7])>=2 && \$ratio>=0.5)){print}' $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt\n";
	print OUT "perl -lane \'if(/lclip_cnt/){print;next} my \@a=split /,/, \$F[8]; my \$ratio=\$a[2]; if(\$F[0] eq \"with_both_clip_disc\" && \$ratio>=0.5){print}' $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/filter_the_bk_with_only_disc_reads.pl $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt2\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/pair_the_end_points.pl $sample/L1/candidate_list_from_disc.txt $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt2 > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/check_ref_repeat_for_somatic_L1.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/ref/repeat_masker/hg38_repeat_masker.tsv.sort $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.add_ref_repeat\n";
#@	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/add_more_candidate_from_vcf.pl $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.add_ref_repeat $sample/L1/$sample\_LINE1.vcf > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine2\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/add_more_candidate_from_vcf.pl $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.add_ref_repeat ../no_threshold_in_clip_and_disc_stp/$sample/L1/$sample\_LINE1.vcf > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine2\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/get_somatic_event_absent_in_all_samples_2nd_round.pl $all_non_ref_list $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine2 $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.absent_from_all_samples.v2_tmp\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/get_somatic_event_absent_in_all_samples_3rd_round.pl $zx_germline_list $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.absent_from_all_samples.v2_tmp $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.absent_from_all_samples.v2\n";
    
	print OUT "\n#add LINE1_cns_pos, and calcuate ins length and tsd length\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/check_polyA_and_cns_pos_score_cigar_from_sam.pl $sample/L1/tmp/$sample.bam.clipped.sam $sample/L1/candidate_list_from_clip.txt_tmp > $sample/L1/candidate_list_from_clip.txt_tmp.add_cns_pos\n";
    print OUT "perl /home/ch237385/bin/BCH/method_development/add_partner_for_single.pl $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.absent_from_all_samples.v2 $sample/L1/candidate_list_from_clip.txt_tmp > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.add_parter_side_for_single\n";
    print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/add_L1_size_tsd_len.pl $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.add_parter_side_for_single $sample/L1/candidate_list_from_clip.txt_tmp.add_cns_pos $sample > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.add_parter_side_add_cns_pos_ins_tsd_size\n";
    print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/label.pl $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt2 $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio.flt.refine.somatic.add_parter_side_add_cns_pos_ins_tsd_size > $sample/L1/ORFeus_insertion.sub-types.txt.refine.somatic.absent_from_all_samples.add_cns_pos.v2_label\n";
	close OUT;
}
close BL;
