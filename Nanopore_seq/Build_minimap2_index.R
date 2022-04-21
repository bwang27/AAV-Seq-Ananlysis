
#setwd("/drive2/wli/analyze/dep_seq/othProject/Nanopore_seq/target_index/")
setwd("/home/Research_Archive/ProcessedData/Nanopore_seq/target_index/")
source("/drive2/wli/analyze/Pipeline/SharedCodes/Rfunc.inc.R")

project_name="2020-06-Research7_CMC5_AAVs_Seq"
project_name="2020-09-Res7-CMC5_AAVSeq"
project_name="2020-10-05_AAVSeq"
project_name="2020-11-03_AAV_SEQ"
project_name="AADC_BioReliance_miSeq"
project_name="2021-01-22_AAV_SEQ"
project_name="2021-03-09_AAV_SEQ_POOL2"
project_name="2021-04-06_AAV_SEQ_MK1C-3CMCTox"
project_name="2021-04-02_AAV_SEQ"
project_name="2021-05-03_AAV_SEQ_CMC_3AS_Justin"

summary_folder="/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/"
plasmid_seq_folder=paste0(summary_folder, "VirusGenomeSequences/")
hg_fa_dir="/drive2/wli/analyze/Pipeline/ReferenceDB/ucsc/genomes/hg19"
Ecoli_fa="Ecoli/GCF_000005845.2_ASM584v2_genomic.fna"
Ad5_fa=paste0(plasmid_seq_folder, "Ad5/Human_Ad5_first5k.fa") 
analysis_folder="/drive2/wli/analyze/dep_seq/othProject/Nanopore_seq/"
sample_info_f=paste0(summary_folder,project_name,"/",project_name,"-Samples.xlsx")
othGOI_fa_dir=paste0(analysis_folder,project_name,"/othGOI_fa/")
out_geneFeat_dir=paste0(analysis_folder,project_name,"/geneFeature/")
E1AB_feature_f=paste0(summary_folder,"/VirusGenomeSequences/Ad5/Human_Ad5_first5k.GOI_info.txt")

library("openxlsx")

sample_info_d=read.xlsx(sample_info_f)
index_d=sample_info_d[,c("pHelper","pRepCap","GOI","GOI_folder")]; dim(index_d)
index_d=index_d[!duplicated(index_d), ]; dim(index_d)

index_d$index_folder=apply(index_d[,1:3],1,paste,collapse=".")
index_d$final_fa_f=paste0(index_d$index_folder,"/",index_d$index_folder,".fa")
index_d$final_mmi_f=paste0(index_d$index_folder,"/",index_d$index_folder,".mmi")
index_d$final_geneFeature_f=paste0(index_d$index_folder,"/",index_d$index_folder,".feature.txt")

plasmid_fa=paste0( plasmid_seq_folder,index_d$GOI_folder,"/",index_d$GOI,".fa ",
	ifelse(is.na(index_d$pHelper)," ", paste0(plasmid_seq_folder,"Helper/", index_d$pHelper, ".fa ")),
	ifelse(is.na(index_d$pRepCap)," ", paste0(plasmid_seq_folder,"RepCap/", index_d$pRepCap, ".fa"))
)
hg_chr_fas=paste0(hg_fa_dir,"/chr[1-9XYM].fa ", hg_fa_dir, "/chr[1-9][0-9].fa " )
index_d$cat_fa_cmd=paste("cat", plasmid_fa, Ad5_fa, Ecoli_fa, hg_chr_fas, ">", index_d$final_fa_f)
index_d$minimap2_cmd=paste("minimap2 -t 8 -d ",index_d$final_mmi_f, index_d$final_fa_f)

##combine gene feature table for all plasmids or selected genes (eg. E1A, E1B)
geneFeature_fs=paste0( plasmid_seq_folder,index_d$GOI_folder,"/",index_d$GOI,".GOI_info.txt ",
	ifelse(is.na(index_d$pHelper)," ", paste0(plasmid_seq_folder,"Helper/", index_d$pHelper, ".GOI_info.txt ")),
	ifelse(is.na(index_d$pRepCap)," ", paste0(plasmid_seq_folder,"RepCap/", index_d$pRepCap, ".GOI_info.txt"))
)
index_d$cat_geneFeature_cmd=paste("cat", geneFeature_fs, E1AB_feature_f, " | perl -ne 'print if ++$k{$_}==1' >", index_d$final_geneFeature_f )

write.table(index_d, file=paste0(sample_info_f,".make_indexes_commands.txt"), col.names=T, row.names=F, sep="\t", quote=F)

forceRewrite=F
for(i in 1:nrow(index_d)){
	if(!file.exists(index_d$final_mmi_f[i]) | forceRewrite){
		print(paste0("making index file: ", index_d$final_mmi_f[i]))
		mkdir_if_not_exist(index_d$final_mmi_f[i])
		system( index_d$cat_fa_cmd[i] )
		system( index_d$minimap2_cmd[i] )
		system(paste("samtools faidx ", index_d$final_fa_f[i]))
		system(paste0("cut -f1,2 ",index_d$final_fa_f[i],".fai > ",index_d$final_fa_f[i],".sizes"))
	}else{
		print(paste0("index file already exist: ", index_d$final_mmi_f[i]))
	}
}

##create a .fa file for all GOI sequences (used for the second step of mapping to remove reads that were wrongly demultiplexed)
all_GOI_fs=unique(paste0( plasmid_seq_folder,index_d$GOI_folder,"/",index_d$GOI,".GOI.fa ")); all_GOI_fs
mkdir_if_not_exist(othGOI_fa_dir)
for(i in 1:nrow(index_d)){
	all_othGOI_fs=setdiff(all_GOI_fs, paste0( plasmid_seq_folder,index_d$GOI_folder,"/",index_d$GOI,".GOI.fa ")[i])
	if(length(all_othGOI_fs)>0){
		out_all_othGOI_fa_f=paste0(othGOI_fa_dir, index_d$index_folder[i],".othGOI.fa")
		cat_othGOI_cmd=paste("cat ", paste(all_othGOI_fs,collapse=" "), " >",out_all_othGOI_fa_f,  sep="")
		print(cat_othGOI_cmd)
		system(cat_othGOI_cmd)
	}
}

##combine gene feature table for all plasmids or selected genes (eg. E1A, E1B)
for(i in 1:nrow(index_d)){
	system( index_d$cat_geneFeature_cmd[i] )
}

