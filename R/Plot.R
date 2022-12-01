library(data.table)

source("R/load_data.R")

##

density.max <- function(x, ... ){
  # calculates the peak position in a density estimate
  if(length(na.omit(x)) < 2){
    return(-Inf)
  }else{
    dw <- density(x, ... )
    dw.max <- dw$x[which.max(dw$y)]
    return(dw.max)
  }
}

##

tableS3 <- fread("TableS3_1-s2.0-S0960982222013811-mmc2.txt")
table66genes <- c("Slc25a5", "Plp2",  "Cox7b", "Trap1a", "Morf4l2", "Ndufb11",  "Renbp",  "Fmr1nb", "Idh3g", "Vbp1", "Ssr4", "Pgk1", "Mageb16", "Rbbp7","Mcts1","Ube2a","Msn","Hsd17b10","1110012L19Rik","Uba1","Ndufa1","Psmd10","Hmgn5","Eif1ax","Mpp1","Pdha1","Praf2","Ldoc1","Igbp1","Mbnl3","Slc6a14","Rnf128","C430049B03Rik","Apoo","Efnb1","Apool","Naa10","Dcaf12l1","Nsdhl","Abcd1","Plac1","Nxf7","Amot","Vma21","Brcc3","Phf6","Klf8","Ubqln2","Ocrl","Gdi1","Slc9a6","Flna","Ccdc22","Xlr3b","Capn6","Zfp275","Slc7a3","Prickle3","Med12","Sept6","Pir","Eda2r","Xist","Fgf13","Srpk3","Tsix")

# Average allelic expression
avg.deng <- deng.allele.melt[expressed == T & chr != "chrY", base::mean(value, na.rm=T, trim = 0.2),by=c("sample","lineage","sex","l1","chrx","maternal")]
avg.borensztein <- borensztein.allele.melt[expressed == T & chr != "chrY", base::mean(value, na.rm=T, trim = 0.2),by=c("sample","genotype","lineage","sex","l1","chrx")]
avg.cheng <- cheng.allele.melt[expressed == T & chr != "chrY", base::mean(value, na.rm=T, trim = 0.2),by=c("sample","lineage", "embryonicday","sex","l1","chrx")]
avg.cheng[chrx == T, c("xi", "xci") := list(V1 == V1[which.min(V1)], abs(0.5 - V1[l1 == "c57"] / sum(V1, na.rm = T))/0.5 ), by="sample"]
avg.lentini <- lentini.allele.melt[expressed == T & chr != "Y", base::mean(value, na.rm=T, trim = 0.2),by=c("sample","day", "sex", "l1", "chrx")]
avg.lentini[chrx == T, c("xi", "xci") := list(V1 == V1[which.min(V1)], abs(0.5 - V1[l1 == "c57"] / sum(V1, na.rm = T))/0.5 ), by="sample"]

# X:A ratios
autoratio.deng <- deng.all.melt[value > 1, median(value[chrx == T], na.rm=T)/median(value[chrx == F],na.rm = T), by=c("sample", "lineage","sex")]
autoratio.borensztein <- borensztein.all.melt[value > 1 & chr != "chrY", median(value[chr == "chrX"], na.rm=T)/median(value[chr != "chrX"], na.rm = T), by=c("genotype", "lineage", "cross", "cell", "sex")]
autoratio.wang <- wang.all.melt[value > 1 & chr %in% paste0("chr", c(1:19,"X")), median(value[chr == "chrX"], na.rm=T)/median(value[chr != "chrX"], na.rm = T), by=c("lineage", "genotype", "sex", "sample")]
autoratio.cheng <- cheng.all.melt[value > 1, median(value[chrx == T], na.rm=T)/median(value[chrx == F],na.rm = T), by=c("sample", "lineage", "embryonicday","sex")]
autoratio.cheng[avg.cheng, xci := xci, on = "sample == sample"]

# Average cluster expression
deng.all.melt[value > 1 & chr %in% paste0("chr", c(1:19,"X")), xa := value/median(value[chr != "chrX"], na.rm = T), by=c("sample")]

cheng.all.melt[value > 1 & chr %in% paste0("chr", c(1:19,"X")), xa := value/median(value[chr != "chrX"], na.rm = T), by=c("sample")]
cheng.all.melt[tableS3, cluster := cluster, on = c("gene == symbol")]
cheng.all.melt[avg.cheng, xci := xci, on = "sample == sample"]

clust.cheng <- cheng.all.melt[!is.na(cluster), median(xa, na.rm=T), by=c("cluster","lineage", "embryonicday", "lineage", "sex", "sample", "xci")]

wang.all.melt[value > 1 & chr %in% paste0("chr", c(1:19,"X")), xa := value/median(value[chr != "chrX"], na.rm = T), by=c("sample")]
wang.all.melt[tableS3, cluster := cluster, on = c("gene == symbol")]

clust.wang <- wang.all.melt[!is.na(cluster), median(xa, na.rm=T), by=c("cluster","lineage", "genotype", "sex", "sample")]

lentini.all.melt[value > 1 & chr %in% c(1:19,"X"), xa := value/median(value[chr != "X"], na.rm = T), by=c("sample")]
lentini.all.melt[tableS3, cluster := cluster, on = c("symbol == symbol")]
lentini.all.melt[avg.lentini, xci := xci, on = "sample == sample"]

clust.lentini <- lentini.all.melt[!is.na(cluster), median(xa, na.rm=T), by=c("cluster", "sex", "day", "sample", "xci")]

# Linear regression
clust.wang[lineage %in% c("X4cell", "X8cell", "E2.5", "E3", "E3.5") & sex == "male", summary(lm(V1~factor(cluster)+lineage+genotype))]
clust.wang[lineage %in% c("X4cell", "X8cell", "E2.5", "E3", "E3.5") & sex == "female" & genotype != "XIST", summary(lm(V1~factor(cluster)+lineage+genotype))]
clust.wang[lineage %in% c("X4cell", "X8cell", "E2.5", "E3", "E3.5") & sex == "female" & genotype != "KO", summary(lm(V1~factor(cluster)+lineage+genotype))]

# Average expression per gene
gn.wang <- wang.all.melt[chr %in% paste0("chr", c(1:19,"X")), mean(value), by=c("gene", "chr", "lineage", "cluster", "genotype")]
gn.wang[, expr.pct := ecdf(V1)(V1), by=c("lineage", "genotype")]

# gn.wang[genotype == "WT" & lineage == "E3.5" & gene %in% table66genes, median(expr.pct)] # expression percentile for the list of 66 genes

# Add annotations
cheng.ratio.melt[avg.cheng, xci := xci, on = "sample == sample"]

cheng.allele.melt[tableS3, cluster := cluster, on = c("gene == symbol")]
cheng.allele.melt[avg.cheng, c("xci", "xi") := list(xci, xi), on = c("sample == sample", "l1 == l1", "chrx == chrx")]

lentini.allele.melt[tableS3, cluster := cluster, on = c("symbol == symbol")]
lentini.allele.melt[avg.lentini, c("xci", "xi") := list(xci, xi), on = c("sample == sample", "l1 == l1", "chrx == chrx")]

##########
## PLOT ##
##########

library(ggplot2)
library(cowplot)
library(ggbeeswarm)

## Autosomal ratios
# Deng
p.xa.deng <-
ggplot(autoratio.deng[!is.na(sex)], aes(x=factor(lineage, levels=c("MIIoocyte", "mid2cell", "4cell", "8cell", "16cell", "earlyblast", "midblast", "lateblast")), y=V1, shape=sex, col="WT")) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=sex)) +
  stat_summary(fun="mean", geom="line", aes(group=sex)) +
  #facet_grid(~sex, scales="free_x", space="free_x", drop=T) +
  labs(x=NULL, y="X:Autosomal ratio", title="Deng et al. 2014 (GSE45719)") +
  coord_cartesian(ylim=c(0.5,1.5)) +
  scale_x_discrete(na.translate=F) +
  scale_color_brewer(palette="Dark2", aesthetics = c("col", "fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Borensztein
p.xa.borensztein <-
ggplot(autoratio.borensztein, aes(x=factor(lineage, levels=c("Oo","2C","4C", "8C", "16C", "32C", "64C")), y=V1, col=genotype, shape=sex)) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(sex,genotype) )) +
  stat_summary(fun="mean", geom="line", aes(group=paste(sex,genotype) )) +
  facet_grid(~sex, scales="free_x", space="free_x") +
  labs(x=NULL, y="X:Autosomal ratio", title="Borensztein et al. 2017 (GSE80810)") +
  coord_cartesian(ylim=c(0.5,1.5)) +
  scale_color_brewer(palette="Dark2", aesthetics = c("col", "fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Wang
p.xa.wang <-
ggplot(autoratio.wang, aes(x=lineage, y=V1, col=genotype, shape=sex)) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(sex,genotype) )) +
  stat_summary(fun="mean", geom="line", aes(group=paste(sex,genotype) )) +
  facet_grid(~sex, scales="free_x", space="free_x") +
  labs(x=NULL, y="X:Autosomal ratio", title="Wang et al. 2016 (GSE101838)") +
  coord_cartesian(ylim=c(0.5,1.5)) +
  scale_color_brewer(palette="Dark2", aesthetics = c("col", "fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Cheng
p.xa.cheng <-
ggplot(autoratio.cheng, aes(x=cut(xci, seq(0,1,0.2)), y=V1, shape=sex, col="WT")) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=sex)) +
  stat_summary(fun="mean", geom="line", aes(group=sex)) +
  facet_grid(~lineage, scales="free_x", space="free_x", drop=T) +
  labs(x=NULL, y="X:Autosomal ratio", title="Chen et al. 2016 (GSE74155) + Cheng et al. 2019 (GSE109071)") +
  coord_cartesian(ylim=c(0.5,1.5)) +
  scale_x_discrete(na.translate=F) +
  scale_color_brewer(palette="Dark2", aesthetics = c("col", "fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))


## Bulk XX, XY, X0
# Werner
p.fc.werner <- 
ggplot(werner.res.melt[!is.na(chrx) & baseMean >= 100],aes(x=log2FoldChange, y=..ndensity.., fill=chrx, col=chrx)) +
  geom_density(alpha=0.33) +
  geom_vline(data=werner.res.melt[chrx==F & baseMean > 100,density.max(log2FoldChange, na.rm=T),by=.id], col="grey", lty=2, aes(xintercept=V1+1)) +
  facet_grid(.id~.) +
  labs(y="Scaled density", x="log2 fold change", title = "Werner et al. 2017 (PRJNA354946)") +
  coord_cartesian(xlim=c(-2,2)) +
  scale_color_brewer(palette="Dark2", name=NULL, labels=c("ChrA","ChrX"), aesthetics = c("col","fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), legend.position = "top", panel.border = element_rect(colour="black"))

## Allelic expression
# Deng
p.ae.deng <-
ggplot(avg.deng[!is.na(sex)], aes(x=factor(lineage, levels=c("MIIoocyte", "mid2cell", "4cell", "8cell", "16cell", "earlyblast", "midblast", "lateblast")), y=V1, col= paste(chrx, l1==tolower(maternal)), shape=sex) ) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(sex,l1==tolower(maternal), chrx) )) +
  stat_summary(fun="mean", geom="line", aes(group=paste(sex,l1==tolower(maternal), chrx) )) +
  facet_grid(~sex, scales="free_x", space="free_x") +
  labs(x=NULL, y="Allelic expression (chrX)", title="Deng et al. 2014 (GSE45719)") +
  coord_cartesian(ylim=c(0,40)) +
  scale_color_brewer(palette="Paired", name="Allele", labels = c("chrA paternal", "chrA maternal", "chrX paternal", "chrX maternal") ) +
  scale_x_discrete(na.translate=F) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Borensztein
p.ae.borensztein <-
ggplot(avg.borensztein[chrx==T & sex =="Female"], aes(x=factor(lineage, levels=c("4C", "8C", "16C", "32C", "64C")), y=V1, col = l1) ) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(genotype, l1))) +
  stat_summary(fun = "mean", geom="line", aes(group=paste(genotype, l1))) +
  facet_grid(~genotype, scales="free_x", space="free_x") +
  labs(x=NULL, y="Allelic expression (chrX)", title="Borensztein et al. 2017 (GSE80810)") +
  coord_cartesian(ylim=c(0,40)) +
  scale_color_brewer(palette="Set1", name="Allele") +
  scale_x_discrete(na.translate=F) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Cheng
p.ae.cheng <-
ggplot(avg.cheng[chrx==T], aes(x= cut(xci, seq(0,1,0.2)), y=V1, col=xi, shape=sex ) ) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(xi, sex) )) +
  stat_summary(fun = "mean", geom="line", aes(group=paste(xi, sex) )) +
  facet_grid(~lineage, scales="free_x", space="free_x") +
  labs(x="rXCI completion (binned)", y="Allelic expression (chrX)", title="Chen et al. 2016 (GSE74155) + Cheng et al. 2019 (GSE109071)") +
  coord_cartesian(ylim=c(0,40)) +
  scale_color_brewer(palette="Set1", name="Allele", labels=c("Xa", "Xi") ) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

## Allelic ratio
# Deng
p.ratio.deng <-
ggplot(deng.ratio.melt[chrx == T & expressed == T & !is.na(sex), median(value, na.rm=T), by=c("sample", "lineage", "sex", "maternal")], aes(x=factor(lineage, levels=c("MIIoocyte", "mid2cell", "4cell", "8cell", "16cell", "earlyblast", "midblast", "lateblast")), y=ifelse(maternal == "C57", V1, 1-V1), shape=sex ) ) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=sex)) +
  stat_summary(fun = "mean", geom="line", aes(group=sex)) +
  geom_hline(yintercept = 0.5, lty=2, col="grey") +
  facet_grid(~sex, scales="free_x", space="free_x") +
  labs(x=NULL, y="Maternal fraction (chrX)", title="Deng et al. 2014 (GSE45719)") +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_discrete(na.translate=F) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Borensztein
p.ratio.borensztein <-
ggplot(borensztein.ratio.melt[chrx == T & sex == "Female" & expressed == T, median(value, na.rm=T), by=c("sample", "lineage", "sex", "genotype")], aes(x=factor(lineage, levels=c("4C", "8C", "16C", "32C", "64C")), y=1-V1) ) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=genotype)) +
  stat_summary(fun = "mean", geom="line", aes(group=genotype)) +
  geom_hline(yintercept = 0.5, lty=2, col="grey") +
  facet_grid(~genotype, scales="free_x", space="free_x") +
  labs(x=NULL, y="Maternal fraction (chrX)", title="Borensztein et al. 2017 (GSE80810)") +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_discrete(na.translate=F) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

## Cluster expression
# Wang
p.cluster.wang <-
ggplot(clust.wang[lineage %in% c("X4cell", "X8cell", "E2.5", "E3", "E3.5")], aes(x=lineage, y=V1, col=genotype, shape=sex )) + 
  stat_summary(fun.data = "mean_cl_boot", geom="ribbon", col=NA, alpha=0.33, aes(group=genotype, fill=genotype) ) + 
  stat_summary(fun = "mean", geom="point") + 
  stat_summary(fun="mean", geom="line", aes(group=genotype )) +
  facet_grid(sex~cluster) +
  labs(x=NULL, y="X:Autosomal ratio", title="Wang et al. 2016 (GSE101838)") +
  coord_cartesian(ylim=c(0,4)) +
  scale_colour_brewer(palette="Dark2", aesthetics = c("col", "fill")) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Cheng
p.cluster.cheng <-
ggplot(clust.cheng[lineage %in% "EPI"], aes(x=cut(xci, seq(0,1,0.2)), y=V1, shape=sex )) + 
  stat_summary(fun.data = "mean_cl_boot", geom="ribbon", col=NA, alpha=0.33, aes(group=sex) ) +
  stat_summary(data = clust.cheng[lineage %in% "EPI" & sex == "M"], fun.data = "mean_cl_boot", show.legend = F, aes(group=sex) ) +
  stat_summary(fun = "mean", geom="point") + 
  stat_summary(fun="mean", geom="line", aes(group=sex )) +
  facet_grid(~cluster) +
  labs(x="rXCI completion (binned)", y="X:Autosomal ratio", title="Chen et al. 2016 (GSE74155) + Cheng et al. 2019 (GSE109071)") +
  coord_cartesian(ylim=c(0,4)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Lentini
p.cluster.lentini <-
ggplot(clust.lentini[!xci %between% c(0.2,0.8)], aes(x=cut(xci, seq(0,1,0.2)), y=V1)) + 
  stat_summary(fun.data = "mean_cl_boot", geom="ribbon", col=NA, alpha=0.33, aes(group=cluster) ) +
  stat_summary(fun="mean", geom="point") + 
  stat_summary(fun="mean", geom="line", aes(group=cluster )) +
  facet_grid(~cluster) +
  labs(x="rXCI completion (binned)", y="X:Autosomal ratio", title="Lentini et al. 2022 (E-MTAB-9324)") +
  coord_cartesian(ylim=c(0,4)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

## Allelic cluster expression
# Cheng
p.clustae.cheng <-
ggplot(cheng.allele.melt[expressed == T & !is.na(cluster) & lineage == "EPI", base::mean(value, na.rm=T, trim=0.2), by=c("sample", "sex", "l1", "xci", "xi", "cluster")], aes(x=cut(xci, seq(0,1,0.2)), y=V1, col=xi, shape=sex )) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(xi, sex) ) ) +
  stat_summary(fun = "mean", geom="line", aes(group=paste(xi, sex) )) +
  facet_grid(~cluster) +
  labs(x="rXCI completion (binned)", y="Allelic expression", title="Chen et al. 2016 (GSE74155) + Cheng et al. 2019 (GSE109071)") +
  scale_color_brewer(palette="Set1", name="Allele", labels=c("Xa", "Xi") ) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Lentini
p.clustae.lentini.rpkm <-
ggplot(lentini.allele.melt[expressed == T & !is.na(cluster) & !xci %between% c(0.2, 0.8), base::mean(value, na.rm=T, trim=0.2), by=c("sample", "sex", "l1", "xci", "xi", "cluster")], aes(x=cut(xci, seq(0,1,0.2)), y=V1, col=xi )) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(xi) ) ) +
  stat_summary(fun = "mean", geom="line", aes(group=paste(xi) )) +
  facet_grid(~cluster) +
  labs(x="rXCI completion (binned)", y="Allelic expression (RPKM)", title="Lentini et al. 2022 (E-MTAB-9324)") +
  scale_color_brewer(palette="Set1", name="Allele", labels=c("Xa", "Xi") ) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

p.clustae.lentini.relumi <-
ggplot(lentini.allele.melt[expressed == T & !is.na(cluster) & !xci %between% c(0.2, 0.8), base::mean(relumi, na.rm=T, trim=0.2), by=c("sample", "sex", "l1", "xci", "xi", "cluster")], aes(x=cut(xci, seq(0,1,0.2)), y=V1, col=xi )) +
  stat_summary(fun.data = "mean_cl_boot", aes(group=paste(xi) ) ) +
  stat_summary(fun = "mean", geom="line", aes(group=paste(xi) )) +
  facet_grid(~cluster) +
  labs(x="rXCI completion (binned)", y="Allelic expression (rel. UMI)", title="Lentini et al. 2022 (E-MTAB-9324)") +
  scale_color_brewer(palette="Set1", name="Allele", labels=c("Xa", "Xi") ) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour="black"), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

## 66 genes expression
# Deng
p.66.deng <-
ggplot(deng.all.melt[gene %in% table66genes & sex == "F" & expressed == T], aes(x=factor(lineage, levels=c("MIIoocyte", "mid2cell", "4cell", "8cell", "16cell", "earlyblast", "midblast", "lateblast")), y=xa)) + 
  stat_summary(fun="median", geom="line", aes(group=1)) +
  stat_summary(fun.data="median_cl_boot") + 
  labs(x="Pre-implantation", y="X:Autosomal ratio", title = "(GSE45719)") +
  coord_cartesian(ylim=c(0,3)) +
  scale_x_discrete(na.translate=F) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

# Cheng
p.66.cheng <-
ggplot(cheng.all.melt[gene %in% table66genes & lineage == "EPI" & sex == "F" & expressed == T], aes(x=cut(xci,seq(0,1,0.2)), y=xa)) + 
  stat_summary(fun="median", geom="line", aes(group=1)) +
  stat_summary(fun.data="median_cl_boot") + 
  labs(x="Post-implantation", y="X:Autosomal ratio", title = "(GSE74155+GSE109071)") +
  coord_cartesian(ylim=c(0,3)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

## Cluster percentiles
# Wang
p.pct.wang <-
ggplot(gn.wang[chr == "chrX" & genotype == "WT" & !is.na(cluster), median(expr.pct), by=c("lineage", "cluster")], aes(x=factor(cluster), y=V1)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0.5, lty=2, col="grey") +
  labs(x="Cluster", y="Expression percentile", title="Wang et al. 2016 (GSE101838)") +
  coord_cartesian(ylim=c(0.25,0.75)) +
  scale_colour_viridis_d(option="turbo") +
  theme_cowplot()
  
## EXPORT

ggsave2("plots/xa_ratio.pdf", width = 18, height = 4,
  plot_grid(nrow = 1, rel_widths = c(1,1.4,1.6, 1),align = "h", axis = "tb",
    p.xa.deng,
    p.xa.borensztein,
    p.xa.wang,
    p.xa.cheng
))

ggsave2("plots/allelic_expression.pdf", width = 10, height = 6,
  plot_grid(nrow = 2, rel_widths = c(1,0.8), rel_heights = c(1, 0.7),align = "hv", axis = "trbl",
    p.ae.deng, p.ae.borensztein,
    p.ratio.deng, p.ratio.borensztein
))

ggsave2("plots/allelic_expression_cheng.pdf", width = 4, height = 3, p.ae.cheng)

ggsave2("plots/bulk_foldchange.pdf",p.fc.werner, width = 3, height = 4)

ggsave2("plots/cluster_xa_ratio.pdf", width = 8, height = 6, 
  plot_grid(nrow=2, rel_heights = c(1,0.8), align="v", axis="trbl",
    p.cluster.wang,
    p.cluster.cheng
))

ggsave2("plots/cluster_allelic_expression.pdf", width = 8, height = 3, 
  p.clustae.cheng
)

ggsave2("plots/cluster_expression_lentini.pdf", width = 4, height = 6,
  plot_grid(nrow=3, align="hv", axis="trbl",
    p.clustae.lentini.rpkm,
    p.clustae.lentini.relumi,
    p.cluster.lentini
))

ggsave2("plots/66gene_xa_ratio.pdf", width = 4, height = 3,
  plot_grid(nrow=1, rel_widths = c(1,0.8), align = "h", axis = "tb", 
    p.66.deng, p.66.cheng
))

ggsave2("plots/cluster_expression_percentile.pdf", width = 3, height = 4, p.pct.wang)
