check_one_dataset<-function(fn,genes) {
  dt<-read.table(fn,header=T,sep="\t",check.names=F,colClass="character")
  good_genes<-intersect(genes,as.character(dt[,1]))
  rlt<-t(dt[dt[,1] %in% good_genes,-c(2:7)])
  rlt2<-cbind(Sample=row.names(rlt),rlt)
  return(rlt2)
}

#setwd("R:/geelegrp/projects/PedDep_DataScience/common/WGS/PedDepAcc/FPKM")


genes<-c("GeneName","MYCN","MYC","PBX1","RUNX1T1","TAL1","TAL2","LMO1","LMO2","LYL1","TLX1","TLX3","NKX2-1",
   "HOXA9","STAG2","BCL11B","NOTCH1","NKX2-5","SPI1","CDK6","CRLF2","MYB","BCL2","EPOR","HLF",
   "TEF","DBP","MDM2","CDK4","DUX4","EGFR","ERBB2","EGFR","CCND1","MET","TERT","LCK","YAP1","CCNE1",
   "CCNE2","ALK","FGFR1","FGFR2","FGFR3","PDGFRA","KIT","NSD2","HOXA3","HOXD3","HOXC11","HOXB3",
   "FOXK2","FOXR2","LIN28B","LCK","CTNNB1","ABL1","CCNA1","CCNA2","CCNB1","CCNB2","CCNB3","CCNC",
   "CCND2","CCND3","CCNF","CCNG1","CCNG2","CCNH","CCNI","CCNI2","CCNJ","CCNJL","CCNK","CCNL1",
   "CCNL2","CCNO","CCNT1","CCNT2","CCNY","CCNYL1","TP53","KMT2D","ARID1A","PTEN","RB1","CREBBP",
   "APC","EP300","NF1","BCORL1","SMARCA4","CDKN2A","MSH6","BRCA2","CTCF","HDAC2","MGA","ARID1B",
   "JAK1","MSH2","MSH3","KDM6A","NF2","NIPBL","BCOR","PBRM1","KMT2A","PTCH1","ATM","SMAD4","STK11",
   "SETD2","TGFBR2","MLH1","ATF7IP","ARID2","B2M","DNMT3A","PIK3R1","PMS2","AMER1","WT1","USP7","NR3C1",
   "POLE","TET2","HLA-A","ETV6","BRCA1","SMAD2","VHL","FAT1","PAX5","ATRX","SMARCB1","SMO","RAD21","GATA1",
   "TSPYL2","CDKN1B","ZBTB7A","DNM2","ASXL2","FBXW7","FAT2","CDKN1A","USP9X","TSC2","EZH2","ID3","SUFU",
   "DICER1","SUZ12","ASXL1","DDX3X","PHF6","GATA3","PPM1D","KEAP1","CASP8","XBP1","POLQ","DROSHA","TGFBR1",
   "MED12","SH2B3","BLM","MAP2K4","IKZF3","TNFAIP3","CDKN2C","RPL5","BTG1","STAT3","ZFP36L2","CXCR4","PMS1",
   "BARD1","BAP1","TBL1XR1","LEF1","NPM1","JAK2","NONO","UTY","PDX1","FLT3","ELF1","PDK1","BCR","IKZF1","CDX2",
   "PHOX2B","KRAS","ACVR1","AKT1","AR","BCL6","BRAF","CARD11","CFH","CHD4","EBF1","ERBB3","FOXP1","GNAS","GNB1",
   "H3F3A","HRAS","IDH1","IDH2","IKBKB","IL7R","JAK3","MAP2K1","MAP3K1","MAP3K4","MAPK1","MAX","MEF2B","MTOR",
   "NFE2L2","NRAS","NT5C2","PIK3CA","PIK3CD","PPP2R1A","PRDM9","PTPN11","RAC1","RHOA","RIT1","RRAS","RRAS2",
   "SF3B1","SMAD3","SPOP","STAT5B","TP63","TRAF3","TRAF7","TSHR","XPO1","ZEB2","ZIC1","ZMYM3","RAF1","MECOM",
   "CTNNB1","OTX2","NUTM1","NTRK3","MNX1","ZNF384","MEF2D","GLIS2","HNRNPUL1","BCL9","CEBPB","ABL2")

genes<-c("ZNF384","MEF2D","GLIS2","HNRNPUL1","BCL9","CEBPB","ABL2")

rlt1<-check_one_dataset("PROPEL_hg38_RNASEQ_all_TPM.txt",genes)

write.table(rlt1,"ZZ_extrated_genes.txt",row.names=F,col.names=F,sep="\t",quote=F)

