#!/usr/bin/perl

use strict;
use warnings;

################################################################################################################################################
# This program is to generate shell to obtain the clip/disc/polyA info only using the output of clip and disc steps of xTea (1st and 2nd step)##
################################################################################################################################################

my $bam_list = shift;

die "perl $0 <bam.list>\n" unless $bam_list;

open BL, $bam_list or die $!;
while(<BL>){
	chomp;
	my @a = split;
	my $sample = $a[0];
	open OUT, ">$sample/L1/get_clip_disc_polyA_Info.sh" or die $!;
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH -c 2\n";
	print OUT "#SBATCH -t 0-8:10\n";
	print OUT "#SBATCH --mem=30G\n";
	print OUT "#SBATCH -p bch-compute\n";
	print OUT "#SBATCH -o $sample\_get_clip_disc_polyA_%j.out\n\n";
    
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/check_polyA_from_sam.pl $sample/L1/tmp/$sample.bam.clipped.sam_partial_polyA.sam $sample/L1/tmp/$sample.bam.clipped.sam.sam $sample/L1/candidate_list_from_clip.txt_tmp > $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/add_disc_info_at_candidate_list_from_clip.txt_tmp.pl $sample/L1/tmp/raw_discordant_reads_tmp0_seplf $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam > $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/add_polyA_info_again.pl $sample/L1_add_polyA_at_clip_stp/candidate_list_from_clip.txt_tmp $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc > $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc.add_xTea_polyA\n";
	print OUT "perl -lane \'next if /merged_in_peak/; if(/lclip_cnt/){print \"type1/2\\t\$_\"; next}; if((\$F[5]+\$F[6])>0 && (\$F[12]+\$F[13])>0){print \"with_both_clip_disc\\t\$_\"}elsif(\$F[5]+\$F[6]>0){print \"with_only_clip\\t\$_\";}elsif(\$F[12]+\$F[13]>0){print \"with_only_disc\\t\$_\";}else{print \"other\\t\$_\"}\' $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc.add_xTea_polyA | sort > $sample/L1/ORFeus_insertion.sub-types.txt\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_ORFeus/cal_realign_cnt_2-types.pl $sample/L1/ORFeus_insertion.sub-types.txt $sample/L1/tmp/$sample.bam.clipped.sam > $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/check_polyA_and_cns_pos_score_cigar_from_sam.pl $sample/L1/tmp/$sample.bam.clipped.sam $sample/L1/candidate_list_from_clip.txt_tmp > $sample/L1/candidate_list_from_clip.txt_tmp.add_cns_pos\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/add_L1_pos_clip_reads_disc_reads.pl $sample/L1/candidate_list_from_clip.txt_tmp.add_cns_pos $sample/L1/tmp/raw_discordant_reads_tmp0_seplf $sample/L1/ORFeus_insertion.sub-types.txt.add_clip_ORFeus_ratio > $sample/L1/clip_disc_polyA_Info.txt\n";
	print OUT "perl -lane 'print \"\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[4]\\t\$F[12]\\t\$F[13]\\t\$F[21]\\t\$F[22]\\t\$F[23]\\t\$F[24]\"' $sample/L1/clip_disc_polyA_Info.txt > $sample/L1/clip_disc_polyA_Info.txt.exclude_L1_info\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/convert_format_for_clip_disc_polyA_Info.pl $sample/L1/clip_disc_polyA_Info.txt.exclude_L1_info > $sample/L1/clip_disc_polyA_Info.txt.exclude_L1_info.ts\n";
	print OUT "/home/ch237385/bin/BCH/bin/msort -k 1 -k n2 $sample/L1/clip_disc_polyA_Info.txt.exclude_L1_info.ts > $sample/L1/clip_disc_polyA_Info.txt.exclude_L1_info.ts.sort\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/extract_the_pseudogene.pl $sample/L1/clip_disc_polyA_Info.txt.exclude_L1_info.ts.sort > $sample/L1/pseudogenes.txt\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/check_ref_repeat.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/ref/repeat_masker/hg38_repeat_masker.tsv.sort $sample/L1/pseudogenes.txt > $sample/L1/pseudogenes_rm_rep.txt\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/flt_L1-del_by_near_clip.pl $sample/L1/candidate_list_from_clip.txt_tmp.add_polyA_from_sam.add_disc.add_xTea_polyA $sample/L1/pseudogenes_rm_rep.txt > $sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip\n";
	print OUT "perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/select_ins_of_left_eq_right.pl $sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip > $sample/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.final\n";
	close OUT;
}
close BL;
