library(data.table)
library(ggplot2)
library(ape)
library(gridExtra)
library(RColorBrewer)
library(reshape)
library(ggrepel)
library(qvalue)
library(scales)

#Filepaths
datapath <- "./functional_gene_content/"
figpath <- "./functional_gene_content/"

Nlist <- c("N01", "N02", "N03", "N04", "N06", "N07", "N08", "N09", "N10", "N11", "N12", "N13", "N14", "N15", "N16")
Flist <- c("F01", "F02", "F03", "F04", "F06", "F07", "F08", "F09", "F10", "F11", "F12", "F13", "F14", "F15", "F16")
samplelist <- c(Nlist, Flist)
phylo_list <- c("api","bapis","bifido","bom","com","firm4","firm5","fper","gilli","lkun","snod")
treatments <- c("Nurses","Foragers")
sdp_list_restricted <- c("bifido_1.1", "bifido_1.2","bifido_1.3","bifido_1.4",
                         "bifido_1.5","bifido_2","firm4_1","firm4_2","firm5_1",
                         "firm5_2","firm5_3","firm5_4","gilli_1","gilli_2",
                         "snod_1","bapis", "fper_1", "com_1")
COGlist <- c("C","E","F","G","H","I","P","Q","D","M","N","O","T","U","V","W","Z","J","K","L")

### Contig analysis

#Loading data
contigs_x_SDP <- fread(paste0(datapath, "contig_x_SDP.txt"),h=T)
contigs_x_SDP[, maxval:=do.call(pmax, .SD[, -c("contig")])]
contigs_x_SDP[, maxSDP:=colnames(contigs_x_SDP[, -c("contig", "maxval")])[apply(contigs_x_SDP[, -c("contig", "maxval")], 1, which.max)]]
contigs_x_SDP[, total:=rowSums(.SD[, -c("contig", "maxval", "maxSDP")])]
unanimous_SDP <- copy( contigs_x_SDP[ !maxval<0.8 | (maxval==total & maxval>0.5),])
nonunanimous_SDP <- copy(contigs_x_SDP[ maxval<=0.5 | (maxval!=total & maxval<0.8),])

nb_contigs_assigned <- copy(unanimous_SDP[, .N, by=maxSDP])
nb_contigs_assigned <- rbind(nb_contigs_assigned,nonunanimous_SDP[,.N,by=maxSDP])

coreness <- fread(paste0(datapath, "core_OGs_list.txt"),h=F)
setnames(coreness, c("OG", "gene", "phylo"))
coreness[, core:=phylo]
orfcov_dt <- data.table()
for (i in samplelist) {
  tmp_dt <- fread(paste0(datapath, i, "_catalogue_cov.txt"), header = TRUE)
  tmp_dt[, Sample := i]
  orfcov_dt <- rbind(orfcov_dt, tmp_dt)
}
setnames(orfcov_dt, c("orf", "ORFstart", "ORFend", "reads", "covbases", "covper", "coverage", "meanbaseq", "meanmapq", "Sample"))
orfcov_dt[, contig:=paste(strsplit(orf, "_")[[1]][-8], collapse = "_"), by=orf]
orfcov_dt[Sample %like% "N[0-1][0-9]", Host:="Nurses"]
orfcov_dt[Sample %like% "F[0-1][0-9]", Host:="Foragers"]
OGs <- fread(paste0(datapath, "all_filt_ffn.blastn.1.filt.SDP.OG"), header = FALSE)
setnames(OGs, c("orf", "contig", "gene", "strain", "phylo", "sdp", "OG"))
orfcov_dt_an <- merge.data.table(orfcov_dt, unanimous_phylo[, c("contig", "maxphylo")], by = "contig", all.x = TRUE)
orfcov_dt_an <- merge.data.table(orfcov_dt_an, unanimous_SDP[, c("contig", "maxSDP")], by = "contig", all.x = TRUE)
orfcov_dt_an <- merge.data.table(orfcov_dt_an, OGs[, c("orf", "OG")], by="orf", all.x = TRUE)
orfcov_dt_an[, phylo:=strsplit(maxSDP,"_")[[1]][1], by="maxSDP"]
orfcov_dt_an <- merge.data.table(orfcov_dt_an, coreness[, .N, by=c("OG","core","phylo")][, c("OG","core","phylo")], by=c("OG","phylo"), all.x = TRUE)
orfcov_dt_an[is.na(core) & !is.na(phylo), core:="no"]
orfcov_dt_an[core==phylo, core:="yes"]

orfcov_clustered_dt <- copy(orfcov_dt_an[!is.na(maxSDP) & !is.na(OG), sum(coverage, na.rm=TRUE), by=c("Sample", "Host", "maxphylo", "maxSDP", "OG", "core")])

#plots

cube_root <- function(x) x^(1/3)
cube_pow <- function(x) x^3
cbrt <- trans_new(name="cbrt", transform=cube_root, inverse = cube_pow)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

scatplotorfs <- list()
volcanopltlist <- list()
allvolcanosigclusters <- c()
orfcov_cast_dt <- data.table()
allsigclust_dt <- data.table()
all_og_normcov_dt <- data.table()
for (i in sdp_list_restricted) {
  tmp_dt <- copy(orfcov_clustered_dt[maxSDP==i, ])
  for (j in unique(tmp_dt[,Sample])) {
    norm_factor <- 10000/tmp_dt[core=="yes" & Sample==j, sum(V1)]
    tmp_dt[Sample==j, normdepth:=V1*norm_factor]
  }
  tmp_dt_c <- tmp_dt[, mean(normdepth, na.rm=TRUE), by=c("Host", "OG","core")]
  tmp_cast_dt <- dcast.data.table(tmp_dt_c, OG + core ~  Host, value.var="V1", fill = 0)
  tmp_cast_dt[, SDP:=i]
  orfcov_cast_dt <- rbind(orfcov_cast_dt, tmp_cast_dt)
  maxcov=max(c(max(tmp_cast_dt$Foragers, na.rm = TRUE), max(tmp_cast_dt$Nurses, na.rm = TRUE)))
  ddev=sd(c(tmp_cast_dt$Foragers, tmp_cast_dt$Nurses), na.rm=TRUE)
  scatplotorfs[[i]] <- ggplot(tmp_cast_dt, aes(x=Nurses, y=Foragers))+
    geom_point(aes(color=core))+
    scale_color_manual(values=c("#000000", "#C77CFF"))+
    scale_x_continuous("Normalized mean coverage in nurses", limits = c(0,maxcov+1))+
    scale_y_continuous("Normalized mean coverage in foragers", limits = c(0,maxcov+1))+
    coord_trans(x=cbrt, y=cbrt)+
    geom_segment(x=0, xend=maxcov+1,y=0,yend=maxcov+1, color="red")+
    theme_set(theme_classic(base_size=16))+
    ggtitle(i)
  a <- dcast.data.table(tmp_dt, OG ~ Sample, value.var = "normdepth", fill=0)
  a[, fc:=0]
  a[, pval:=0]
  for (j in 1:nrow(a)) {
    n <- as.numeric(melt.data.table(a[j, Nlist[Nlist %in% names(a)],with=FALSE], 
                                    measure.vars = Nlist[Nlist %in% names(a)])$value)
    f <- as.numeric(melt.data.table(a[j, Flist[Flist %in% names(a)],with=FALSE], 
                                    measure.vars = Flist[Flist %in% names(a)])$value)
    a[j, pval:=wilcox.test(f,n, paired=TRUE)$p.value]
    a[j, fc:=log2(mean(f)/mean(n))]
  }
  a[, qval:=qvalue(a[, pval])$qvalues]
  a[, sig := "Non-significant"]
  a[, SDP := i]
  a[as.numeric(qval)<0.05 & as.numeric(fc)>1, sig:= "Foragers"]
  a[as.numeric(qval)<0.05 & as.numeric(fc)<(-1), sig:="Nurses"]
  volcanopltlist[[i]]<-ggplot(a, aes(x=as.numeric(fc), y=as.numeric(qval)))+
    geom_point(aes(color=factor(sig, 
                                levels = c("Nurses", "Foragers", "Non-significant"))),
               alpha=0.5)+
    scale_color_manual(name="Enrichment", values=c("Nurses"="#D55E00",
                                                   "Foragers"="#56B4E9",
                                                   "Non-significant"="#000000"))+
    scale_y_continuous(name="q-value", trans=reverselog_trans(10))+
    scale_x_continuous(name="log2 fold change(Foragers/Nurses)")+
    geom_vline(xintercept=0, linetype="longdash", color="grey")+
    geom_vline(xintercept=-1, linetype="dotted", color="black")+
    geom_vline(xintercept=1, linetype="dotted", color="black")+
    geom_hline(yintercept=0.05, linetype="longdash", color="black")+
    theme_set(theme_classic(base_size=16))+
    ggtitle(i)
  allvolcanosigclusters <- c(allvolcanosigclusters, a[sig=="Nurses" | sig=="Foragers", OG])
  allsigclust_dt <- rbind(allsigclust_dt, a[sig=="Nurses" | sig=="Foragers", c("OG", "fc", "pval", "qval", "sig", "SDP")])
  all_og_normcov_dt <- rbind(all_og_normcov_dt, a[, c("OG", "fc", "pval", "qval", "sig", "SDP")])
}
#orfcov_clustered_dt[, .N, by=c("OG", "maxSDP")][, .N, by="OG"][, .N, by="N"]
pdf(paste0(figpath, "OG_db_orfcovOG.pdf"), onefile=TRUE, width=12)
marrangeGrob(scatplotorfs, nrow = 1, ncol = 1)
dev.off()
write(orfcov_cast_dt[core=="yes" & (Nurses==0 | Foragers==0),]$OG, file=paste0(datapath, "zerocovOGs.txt"), ncolumns = 1, append=FALSE)
pdf(paste0(figpath, "OG_db_volcanoplotsOG.pdf"), onefile=TRUE, width=12)
marrangeGrob(volcanopltlist, nrow = 1, ncol = 1)
dev.off()


annot_all_OGs <- fread(paste0(datapath, "all_ogs_eggnog_annot.emapper.annotations"),h=F,sep="\t")
annot_all_OGs[, V1:=gsub(".fa","",V1)]
annot_all_OGs[, OG:=strsplit(V1,":")[[1]][2], by=V1]
tmp <- copy(annot_all_OGs[, .N, by=c("OG","V21","V22")])
tmp[, idx:=frankv(copy(.SD), cols=c("N"), order=-1, na.last = TRUE, ties.method="random"), by=OG] #Majority vote on the annotation
annot_all_OGs_s <- copy(tmp[idx==1, c("OG","V21","V22")])
setnames(annot_all_OGs_s, c("OG","COG","Annotation"))
annot_all_OGs_we <- merge.data.table(annot_all_OGs_s, all_og_normcov_dt, by=c("OG"), all = TRUE)
annot_all_OGs_we[is.na(COG), COG:=""]
rm(tmp)

annot_sig_OGs_we <- merge.data.table(allsigclust_dt, annot_all_OGs_we[, c("OG","SDP","sig","COG","Annotation")], all.x = TRUE, by=c("OG","SDP","sig"))
write.table(annot_sig_OGs_we, file=paste0(datapath, "allsignifOGs.txt"),append = FALSE, sep = "\t",eol="\n",quote=FALSE)

fishtest_reslist <- data.table()
for (i in sdp_list_restricted){
  for (j in COGlist){
    tmp <-  matrix(c(annot_all_OGs_we[SDP==i & COG %like% j, .N], annot_sig_OGs_we[SDP==i & COG %like% j & sig=="Foragers", .N],
                     annot_all_OGs_we[SDP==i, .N], annot_sig_OGs_we[SDP==i & sig=="Foragers", .N]), nrow = 2,dimnames = list(c("all", "sig"), c(j, "other")))
    fishtest_reslist<-rbind(fishtest_reslist, data.table(i,j,fisher.test(tmp)$p.value,fisher.test(tmp)$estimate))
  }
}
setnames(fishtest_reslist, c("SDP","COG","pvalue","oddsratio"))
fishtest_reslist[pvalue>1, pvalue:=1]
fishtest_reslist[, qval:=qvalue(fishtest_reslist[,pvalue])$qvalues]

#COG categories for significant OGs.

pdf(paste0(figpath, "sigOG_nb_by_SDP.pdf"), onefile=TRUE, width=12)
ggplot(annot_sig_OGs_we[, .N, by=c("OG", "SDP", "sig")][, .N, by=c("SDP","sig")][SDP %in% sdp_list_restricted,], 
       aes(x=SDP, y=N, fill=factor(sig, levels = c("Nurses", "Foragers"))))+
  geom_col(position=position_dodge(preserve = "single"))+
  scale_fill_manual(name="Enrichment", values=c("Nurses"="#D55E00",
                                                "Foragers"="#56B4E9"))+
  scale_y_continuous(name="Number of significant OGs")+
  theme_gray(base_size=22)+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=16),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust = 0.5))
dev.off()


tmp <- annot_sig_OGs_we[, .N, by=c("OG", "SDP", "COG")][, .N, by=c("SDP","COG")][SDP %in% sdp_list_restricted,][, lapply(.SD, function(x) unlist(tstrsplit(x, "", fixed=TRUE))), .SDcols="COG", by =c("SDP","COG","N")][, -2]
tmp[, prop := N / sum(N), by=SDP]
tmp[COG!="S" & COG!="NA", propnoS := N / sum(N), by=SDP]
tmp2 <- annot_sig_OGs_we[SDP %in% sdp_list_restricted, .N, by=c("OG", "COG")][, .N, by=c("COG")][, lapply(.SD, function(x) unlist(tstrsplit(x, "", fixed=TRUE))), .SDcols="COG", by =c("COG","N")][,-1][, sum(N), by=c("COG")]
setnames(tmp2, c("COG","N"))
tmp2[, prop := N / sum(N)]
tmp2[COG!="S" & COG!="NA", propnoS := N / sum(N)]
tmp2[, SDP := "all"]
tmp3 <- rbind(tmp,tmp2)
rm(tmp)
rm(tmp2)

pdf(paste0(figpath, "sigOG_prop_COG.pdf"), onefile=TRUE, width=12)
ggplot(tmp3[COG!="S" & COG!="NA",], aes(x=SDP, y=propnoS, fill=factor(COG, levels=c("C","E","F","G","H","I","P","Q","D","M","N","O","T","U","V","W","Z","J","K","L"))))+
  geom_col()+
  scale_y_continuous(name="Number of significant OGs")+
  scale_fill_manual(name="COG", values=c("C"="#f5eb94","E"="#dea647","F"="#d67f3a","G"="#8b1022",
                                         "H"="#edf6e3","I"="#bbdea4","P"="#81b867","Q"="#457939",
                                         "D"="#3F007D","M"="#54278F","N"="#6A51A3","O"="#807DBA",
                                         "T"="#9E9AC8","U"="#BCBDDC","V"="#DADAEB","W"="#EFEDF5","Z"="#FCFBFD",
                                         "J"="#bfbfbf","K"="#838383","L"="#414141"))+
  theme_set(theme_classic(base_size=22))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=16),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust = 0.5))
dev.off()
