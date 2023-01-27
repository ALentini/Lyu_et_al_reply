library(data.table)
library(magrittr)

### set up functions

fpkm.scale <- function(fpkm, scale.factor){
  fpkm * (scale.factor / sum(fpkm, na.rm=T))
}

scale.allele <- function(x, ratio, labs = c("c57", "cast")){
  idx.row <- intersect(row.names(x), row.names(ratio))
  idx.col <- intersect(colnames(x), colnames(ratio))
  
  out <- list(
    as.matrix(x[idx.row, idx.col]) * as.matrix(ratio[idx.row, idx.col]),
    as.matrix(x[idx.row, idx.col]) * (1 - as.matrix(ratio[idx.row, idx.col]))
  )  
  names(out) <- labs
  return(out)
}

add.info <- function(x, rowdata, coldata, rowjoin = "Var1==name2", coljoin = "Var2==LibraryName"){
  require(data.table)
  require(magrittr)
  
  x %<>% as.data.table()
  rowdata %<>% as.data.table()
  coldata %<>% as.data.table()
  
  x[rowdata, chr := chrom, on = rowjoin]
  x[chr != "chrY" | chr != "Y",chrx := chr == "chrX" | chr == "X"]
  x %<>% .[coldata, on = coljoin]
  
  idx <- match(c(gsub("==.*","",rowjoin),gsub("==.*","",coljoin)),colnames(x))
  
  colnames(x)[idx] <- c("gene", "sample")
  colnames(x) %<>% tolower()
  return(x[!is.na(gene)])
}

loom2matrix <- function(x){
  require(loomR)
  lfile <- connect(x)
  out <- t(lfile$matrix[,])
  colnames(out) <- lfile$col.attrs$CellID[]
  row.names(out) <- lfile$row.attrs$Gene[]
  lfile$close_all()
  return(out)
}

is.outlier <- function(x, nmad=3, ... ){
  x < (median(x, na.rm=T-nmad*mad(x, na.rm=T, ... )))
}

##

## Deng et al. 2014
# load SS1 data
deng.dge.ss1 <- list(
  "rpkm" = as.matrix(read.delim("Deng/ss1_rpkm_all.tsv", row.names = 1, check.names = F)), # check.names=F fixes names starting with number, e.g. 16cell
  "counts.c57" = as.matrix(read.delim("Deng/ss1_counts_c57.tsv", row.names = 1,check.names = F)), 
  "counts.cast" = as.matrix(read.delim("Deng/ss1_counts_cast.tsv", row.names = 1, check.names = F)),
  "meta" = fread("Deng/ss1_metadata.tsv"),
  "refseq" = fread("Deng/ss1_refseq_gene_annotation.tsv")
)
deng.dge.ss1$ratio <- with(deng.dge.ss1, counts.c57 / (counts.c57+counts.cast))
deng.dge.ss1$meta[, c("Lineage", "Maternal") := list(factor(Lineage, levels=unique(Lineage)), gsub(" .*","",GeneticBackground) )]
deng.dge.ss2$rpkm.fix <- apply(deng.dge.ss2$rpkm,2,fpkm.scale, scale.factor = 45e4) # rescale uneven rpkm depth
deng.dge.ss2$allele <- with(deng.dge.ss2, scale.allele(rpkm.fix, ratio))

# load SS2 data
deng.dge.ss2 <- list(
  "rpkm" = as.matrix(read.delim("Deng/ss2_rpkm_all.tsv", row.names = 1, check.names = F)),
  "counts.c57" = as.matrix(read.delim("Deng/ss2_counts_c57.tsv", row.names = 1,check.names = F)),
  "counts.cast" = as.matrix(read.delim("Deng/ss2_counts_cast.tsv", row.names = 1, check.names = F)),
  "meta" = fread("Deng/ss2_metadata.tsv"),
  "refseq" = fread("Deng/ss2_refseq_gene_annotation.tsv")
)
deng.dge.ss2$ratio <- with(deng.dge.ss2, counts.c57 / (counts.c57+counts.cast))
deng.dge.ss2$meta[, c("Lineage", "Maternal") := list(factor(Lineage, levels=unique(Lineage)), gsub(" .*","",GeneticBackground) )]
deng.dge.ss2$tpm <- apply(deng.dge.ss2$rpkm,2,fpkm2tpm) # convert rpkm to tpm
deng.dge.ss2$allele <- with(deng.dge.ss2, scale.allele(rpkm, ratio))

# melt data
deng.all.melt <- rbindlist(list(
  add.info(as.data.table(reshape2::melt(deng.dge.ss1$rpkm.fix)),deng.dge.ss1$refseq, deng.dge.ss1$meta), 
  add.info(as.data.table(reshape2::melt(deng.dge.ss2$rpkm.fix)),deng.dge.ss2$refseq, deng.dge.ss2$meta)
), fill = T)
deng.all.melt[,expressed := mean(value,na.rm = T)>0, by = c("gene","lineage")]
deng.all.melt[lineage == "MIIoocyte", sex := "F"]

deng.allele.melt <- rbindlist(list(
  add.info(as.data.table(reshape2::melt(deng.dge.ss1$allele)),deng.dge.ss1$refseq, deng.dge.ss1$meta), 
  add.info(as.data.table(reshape2::melt(deng.dge.ss2$allele)),deng.dge.ss2$refseq, deng.dge.ss2$meta)
), fill = T)
deng.allele.melt[deng.all.melt, expressed := expressed, on = c("gene", "sample")]
deng.allele.melt[lineage == "MIIoocyte", sex := "F"]

deng.ratio.melt <- rbindlist(list(
  add.info(as.data.table(reshape2::melt(deng.dge.ss1$ratio)),deng.dge.ss1$refseq, deng.dge.ss1$meta), 
  add.info(as.data.table(reshape2::melt(deng.dge.ss2$ratio)),deng.dge.ss2$refseq, deng.dge.ss2$meta)
), fill = T)
deng.ratio.melt[deng.all.melt, expressed := expressed, on = c("gene", "sample")]
deng.ratio.melt[lineage == "MIIoocyte", sex := "F"]

## Borensztein et al. 2017
# load data
borensztein.dge <- list(
  "rprt" = fread("Borensztein/GSE80810_RPRT.txt"), # first two columns are "gene" and "chromosome" so load as data.table
  "ratio" = fread("Borensztein/GSE80810_AllelicRatio.txt"),
  "meta" = as.data.table(t(fread("Borensztein/GSE80810_series_matrix.txt.gz", skip = 32, nrows = 13, header=F)))[-1] # format series matrix
)
borensztein.dge$meta[,sample_id := gsub(".*\\[|\\]","",V1)] %>% .[!grepl("KO", sample_id), sample_id := paste0("WT_",sample_id)]
borensztein.dge$meta[,sex := gsub(".* ","", V13)]
borensztein.dge$allele <- with(borensztein.dge, scale.allele(data.frame(rprt[,-2], row.names = 1), data.frame(ratio[,-2], row.names = 1), labs=c("paternal", "maternal")))

# melt data
borensztein.all.melt <- as.data.table(melt(borensztein.dge$rprt, id.vars=c("Genes","Chromosomes") ))
colnames(borensztein.all.melt)[1:3] <- c("gene", "chr", "sample")
borensztein.all.melt[chr != "chrY",chrx := chr == "chrX"]
borensztein.all.melt[borensztein.dge$meta, sex := sex, on = "sample == sample_id"]
borensztein.all.melt[,c("genotype", "lineage","cross","cell") := tstrsplit(sample,"_")]
borensztein.all.melt[,expressed := mean(value,na.rm = T)>0, by = c("gene","lineage")]

borensztein.allele.melt <- as.data.table(reshape2::melt(borensztein.dge$allele))
colnames(borensztein.allele.melt)[1:4] <- c("gene", "sample", "value", "l1")
borensztein.allele.melt[borensztein.all.melt, c("chr", "sex", "genotype", "lineage", "cross", "cell", "expressed", "chrx") := list(chr, sex, genotype, lineage, cross, cell, expressed, chrx), on = c("gene", "sample")]

borensztein.ratio.melt <- as.data.table(reshape2::melt(borensztein.dge$ratio))
colnames(borensztein.ratio.melt)[1:3] <- c("gene", "chr", "sample")
borensztein.ratio.melt[borensztein.all.melt, c("sex", "genotype", "lineage", "cross", "cell", "expressed", "chrx") := list(sex, genotype, lineage, cross, cell, expressed, chrx), on = c("gene", "sample")]


## Wang et al. 2016
# load and melt data
wang.all.melt <- as.data.table(reshape2::melt(lapply(list.files("Wang", ".txt$", full.names = T), fread)))
colnames(wang.all.melt)[1:2] <- c("gene", "sample")
wang.all.melt[unique(fread("Wang/hgTables")), chr := `#chrom`, on = "gene == name2"]
wang.all.melt[, c("lineage", "genotype", "sex", "replicate") := tstrsplit(sample, "_")]
wang.all.melt[, genotype := toupper(genotype)]
wang.all.melt[, lineage := factor(lineage, levels=c("X4cell", "X8cell", "E2.5", "E3", "E3.5", "E4", "E4.5", "troph"))]
wang.all.melt[, tpm := fpkm2tpm(value), by="sample"]

## Werner et al. 2017
#load data
werner.fls <- list.files(c("werner"),"quant.sf",recursive = T,full.names = T)
names(werner.fls) <- sapply(strsplit(werner.fls,"/"),"[[",2)

library(tximport)
werner.dge <- list(
  "txi" = tximport(werner.fls,type="salmon",tx2gene = fread("gencode.vM22.metadata.MGI.gz", header = F)[,1:2,with=F]),
  "meta" = fread("werner/werner_meta.tsv")
)
werner.dge$meta[,genotype := relevel(factor(genotype),ref="XX")]

# differential expression
library(DESeq2)
werner.dge$dds <- with(werner.dge, DESeqDataSetFromTximport(txi, meta, ~genotype))
werner.dge$dds %<>%  DESeq(test="LRT", reduced = ~1)

werner.res.melt <- rbindlist(idcol=T, list(
  "XXvXO" = as.data.table(as.data.frame(results(werner.dge$dds, contrast = c("genotype", "XX", "XO"))),keep.rownames = T),
  "XXvXY" = as.data.table(as.data.frame(results(werner.dge$dds, contrast = c("genotype", "XX", "XY"))),keep.rownames = T)
))

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
werner.res.melt[, chr := with(getBM(attributes = c("mgi_symbol", "chromosome_name"), filters = "mgi_symbol", values = rn, mart = mart), paste0("chr", chromosome_name)[match(rn, mgi_symbol)] ) ]
werner.res.melt[ !chr %in% c("chrY", "chrMT"), chrx := chr == "chrX"]

## Chen et al. 2016 + Cheng et al. 2019
# load data
cheng.dge <- list(
  "rpkm" = as.matrix(read.delim("ChenCheng/rpkm_all.tsv", row.names = 1, check.names = F)),
  "counts.c57" = as.matrix(read.delim("ChenCheng/counts_c57.tsv", row.names = 1,check.names = F)), 
  "counts.cast" = as.matrix(read.delim("ChenCheng/counts_cast.tsv", row.names = 1, check.names = F)),
  "meta" = fread("ChenCheng/metadata.tsv"),
  "refseq" = fread("ChenCheng/refseq_gene_annotation.tsv")
)
cheng.dge$ratio <- with(cheng.dge, counts.c57 / (counts.c57+counts.cast))
cheng.dge$meta[, "Maternal" := gsub(" .*","",GeneticBackground)]
cheng.dge$allele <- with(cheng.dge, scale.allele(rpkm, ratio))

# melt data
cheng.all.melt <- add.info(as.data.table(reshape2::melt(cheng.dge$rpkm)),cheng.dge$refseq, cheng.dge$meta)
cheng.all.melt[,expressed := mean(value,na.rm = T)>0, by = c("gene","lineage", "embryonicday")]

cheng.allele.melt <- add.info(as.data.table(reshape2::melt(cheng.dge$allele)),cheng.dge$refseq, cheng.dge$meta)
cheng.allele.melt[cheng.all.melt, expressed := expressed, on = c("gene", "sample")]

cheng.ratio.melt <- add.info(as.data.table(reshape2::melt(cheng.dge$ratio)),cheng.dge$refseq, cheng.dge$meta)
cheng.ratio.melt[cheng.all.melt, expressed := expressed, on = c("gene", "sample")]

## Lentini et al. 2022
# load data
lentini.dge <- list(
  "rpkm" = loom2matrix("Lentini/mESC/zUMIs_output/expression/mESC_diff_clone5.rpkm.exon.all.loom"),
  "umi" = loom2matrix("Lentini/mESC/zUMIs_output/expression/mESC_diff_clone5.umicount.exon.all.loom"),
  "counts.c57" = read.delim("Lentini/mESC/zUMIs_output/allelic/mESC_diff_clone5.BL6_reads.txt", row.names = 1),
  "counts.cast" = read.delim("Lentini/mESC/zUMIs_output/allelic/mESC_diff_clone5.CAST_reads.txt", row.names = 1),
  "meta" = fread("Lentini/mESC/metadata_ss3_clone5diff.tsv"),
  "ensembl" = as.data.table(rtracklayer::import("Mus_musculus.GRCm38.97.chr.gtf.gz", feature.type="gene"))
)
colnames(lentini.dge$ensembl)[1] <- "chrom"
lentini.dge$relumi <- apply(lentini.dge$umi, 2, function(x) (x / sum(x, na.rm=T)) * median(colSums(lentini.dge$umi[,!is.outlier(colSums(lentini.dge$rpkm>0))])) )
lentini.dge$ratio <- with(lentini.dge, counts.c57 / (counts.c57+counts.cast))
lentini.dge$allele$rpkm <- with(lentini.dge, scale.allele(rpkm, ratio))
lentini.dge$allele$relumi <- with(lentini.dge, scale.allele(relumi, ratio))

# melt data
lentini.all.melt <- add.info(as.data.table(reshape2::melt(lentini.dge$rpkm[,!is.outlier(colSums(lentini.dge$rpkm>0))])), lentini.dge$ensembl, lentini.dge$meta, rowjoin = "Var1==gene_id", coljoin = "Var2==sample_bc")
lentini.all.melt[lentini.dge$ensembl, symbol := gene_name, on = "gene==gene_id"]
lentini.all.melt[as.data.table(reshape2::melt(lentini.dge$relumi))[, list(Var1, Var2, relumi=value)], relumi := relumi, on = c("gene==Var1", "sample==Var2")]
lentini.all.melt[,expressed := mean(value,na.rm = T)>0, by = c("gene", "day")]

lentini.allele.melt <- add.info(as.data.table(reshape2::melt(lentini.dge$allele$rpkm)), lentini.dge$ensembl, lentini.dge$meta, rowjoin = "Var1==gene_id", coljoin = "Var2==sample_bc")
lentini.allele.melt %<>% .[sample %in% lentini.all.melt[,unique(sample)]]
lentini.allele.melt[lentini.dge$ensembl, symbol := gene_name, on = "gene==gene_id"]
lentini.allele.melt[as.data.table(reshape2::melt(lentini.dge$allele$relumi))[, list(Var1, Var2, relumi=value, L1)], relumi := relumi, on = c("gene==Var1", "sample==Var2", "l1==L1")]
lentini.allele.melt[lentini.all.melt, expressed := expressed, on = c("gene", "sample")]
