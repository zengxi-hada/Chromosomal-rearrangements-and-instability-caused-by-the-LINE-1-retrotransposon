#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-8:10
#SBATCH --mem=200G
#SBATCH -p bch-compute
#SBATCH -o L1_CL_DoxOrf_DMSO_E12_get_somatic_%j.out

perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/flt_sideRetro_set_v2_for_DMSO.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/project/cell_clonal_genome_instability/pseudogene/DNA-seq/merge/L1_CL_DoxOrf_DMSO_E12/out.vcf /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/L1_CL_DoxOrf_DMSO_E12/L1/pseudogenes_rm_rep.txt.flt_by_near_clip > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/L1_CL_DoxOrf_DMSO_E12/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.intesect_with_sideretro
perl /home/ch237385/bin/BCH/method_development/extend_to_pseudogene/get_somatic_pseudogene.pl /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/all_Dox_exclude_DMSO.list /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/L1_CL_DoxOrf_DMSO_E12/L1/pseudogenes_rm_rep.txt.flt_by_near_clip.intesect_with_sideretro > /lab-share/Gene-Lee-e2/Public/home/xizeng/project/xTea/pseudogene_insertion/no_threshold_in_clip_and_disc_stp/L1_CL_DoxOrf_DMSO_E12/L1/pseudogenes_somatic_v2.txt
