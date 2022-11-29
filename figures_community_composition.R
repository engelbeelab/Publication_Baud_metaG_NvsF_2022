library(data.table)
library(ggplot2)
library(ape)
library(gridExtra)
library(RColorBrewer)
library(reshape)
library(ggrepel)
library(qvalue)
library(vegan)

#Filepaths
datapath_raw <- "./"
datapath <- "./community_composition/"
figpath <- "./community_composition/"
#Set the orders for the plots.
sample_order <- c('N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16',
                  'F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16')
phylo_list <- c("api","bapis","bifido","bom","com","firm4","firm5","fper","gilli","lkun","snod")
phylo_order <- c("firm5","bifido","firm4","gilli","snod","api","fper","bapis","com","lkun","bom")
phylo_order_rv <- rev(phylo_order)
treatments <- c("Nurses","Foragers")
#Set the colors.
set1 <- brewer.pal(9,"Set1")
phylo_col <- c(set1, "light grey", "#DFFF00")
phylo_col_rv <- rev(phylo_col)
treat_colors <- c("#D55E00","#56B4E9","#F0E442")


### Gut Weight comparisons

samp_dt <- fread(paste0(datapath_raw, "SamplingGilles2019.csv"))
samp_dt <- samp_dt[, c("Sample", "Sample_type", "Location", "Number_guts", "Average_gut_mass")]

pdf(paste0(figpath, "gutweightplot.pdf"),onefile=TRUE,width=13)
ggplot(samp_dt,
       aes(x=factor(Sample_type,levels=treatments),
           y = Average_gut_mass,
           fill=factor(Sample_type,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Enrichment", values=c("Nurses"="#D55E00",
                                                "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Sample_type,levels=treatments)),
             shape=22,
             size = 6,
             position=position_dodge2(width=0.6))+
  labs(subtitle=paste0("Wilcoxon sign rank test: \n P = ",
                       wilcox.test(samp_dt[Sample %like% "N", Average_gut_mass],
                                   samp_dt[Sample %like% "F", Average_gut_mass])$p.value),
       fill="Host")+
  ylab("Average gut weight (g)")+
  scale_y_continuous(trans = "log10")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))
dev.off()
  
### Bacterial loads comparisons

#qPCR standard curve values from Kesnerova et al. (2017)
actin_intercept = 35.0942119870867
actin_slope = -3.2699442250388
UV_intercept = 36.5821936471122
UV_slope = -3.35085896083287

filepath <- paste0(datapath_raw, "20190813_AllHives.csv")
CT_dt <- fread(filepath, skip="Well Position")
CT_dt <- fread(filepath, skip="Well Position", nrows = CT_dt[,sum(CT != '')])
CT_dt <- merge(CT_dt[, mean(as.numeric(CT), na.rm=TRUE), keyby=c("Sample Name", "Target Name")], 
               CT_dt[, sd(as.numeric(CT), na.rm=TRUE), keyby=c("Sample Name", "Target Name")])
setnames(CT_dt, c("Sample", "Target", "AverageCT", "SD"))
CT_dt <- CT_dt[ !(Sample %like% "[FN]05" ), ]
CT_dt <- CT_dt[ !(Sample %like% "control" ), ]

CT_dt[ Target %like% "Actin", copies := 10^((AverageCT-actin_intercept)/actin_slope)]
CT_dt[ Target %like% "UV", copies := 10^((AverageCT-UV_intercept)/UV_slope)]

CT_dt[ Sample %like% "N", Host := "Nurses" ]
CT_dt[ Sample %like% "F", Host := "Foragers" ]

sampling_dt <- fread(paste0(datapath, "SamplingGilles2019.csv"))
CT_dt <- merge.data.table(CT_dt, sampling_dt[,c("Sample", "DNA_yield")], by = "Sample", all.x = TRUE)
CT_dt[, DNA := copies*DNA_yield/10]
CT_dt[, Hive:=gsub("[NF]", "", Sample)]
CT_dt[, copiesperngDNA:=copies/10]

pdf(paste0(figpath,"bacloadplot.pdf"),onefile=TRUE,width=8)
ggplot(CT_dt[ Target %like% "UV"],
       aes(x=factor(Host,levels=treatments),
           y=copiesperngDNA,
           fill=factor(Host,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                                "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host,levels=treatments),
                 shape=factor(Host,levels=treatments)), 
             size=6,
             position=position_dodge2(width=0.6)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n p = ", 
                         wilcox.test(CT_dt[ Target %like% "UV" & Sample %like% "N", copiesperngDNA], 
                                     CT_dt[ Target %like% "UV" & Sample %like% "F", copiesperngDNA],
                                     paired=TRUE)$p.value),
       fill="Host",
       shape="Host")+
  ylab("16S copy number per ng of DNA")+
  scale_y_continuous()+
  scale_shape_manual(values=c(21,24))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5), 
        plot.subtitle=element_text(hjust=0.5),
        axis.title.x=element_blank())
dev.off()
actin_copy_bplot <- ggplot(CT_dt[ Target %like% "Actin"], 
                        aes(x=factor(Host,levels=treatments), 
                            y=copies, 
                            fill=factor(Host,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host,levels=treatments)), 
             shape=22,
             size=6)+
  geom_line(aes(group=Hive)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n P = ", 
                         wilcox.test(CT_dt[ Target %like% "Actin" & Sample %like% "N", copies], 
                                     CT_dt[ Target %like% "Actin" & Sample %like% "F", copies],
                                     paired=TRUE)$p.value),
       fill="Sample type")+
  ylab("Actin copy number")+
  scale_y_continuous()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))
uv_DNA_bplot <- ggplot(CT_dt[ Target %like% "UV"], 
                        aes(x=factor(Host,levels=treatments), 
                            y=DNA, 
                            fill=factor(Host,levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host,levels=treatments)), 
             shape=22,
             size=6)+
  geom_line(aes(group=Hive)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n P = ", 
                         wilcox.test(CT_dt[ Target %like% "UV" & Sample %like% "N", DNA], 
                                     CT_dt[ Target %like% "UV" & Sample %like% "F", DNA],
                                     paired=TRUE)$p.value),
       fill="Sample type")+
  ylab("Total estimated 16S gene copies in sample")+
  scale_y_continuous()+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

CT_ratio_dt <- merge.data.table(CT_dt[ Target %like% "Actin",],
                                CT_dt[ Target %like% "UV",],
                                by=c("Sample"))
CT_ratio_dt[, ratio := copies.y / copies.x]
#CT_ratio_dt[, DNAratio := DNA.y / DNA.x]
#CT_ratio_dt <- CT_ratio_dt[, c("Sample", "ratio")]

pdf(paste0(figpath, "bacloadplot.pdf"), onefile=TRUE, width = 8)
ggplot(CT_ratio_dt, 
       aes(x=factor(Host.y, levels=treatments),
           y=ratio,
           fill=factor(Host.y, levels=treatments)))+
  geom_boxplot()+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  geom_point(aes(fill=factor(Host.y, levels=treatments)), 
             shape=22,
             size=6,
             position=position_dodge2(width=0.6)) +
  #geom_line(aes(group = Hive)) +
  labs(subtitle = paste0("Wilcoxon sign rank test: \n P = ",
                         wilcox.test(CT_ratio_dt[Sample %like% "N", ratio],
                                     CT_ratio_dt[Sample %like% "F", ratio])$p.value),
       fill="Sample type")+
  scale_y_continuous(name="16S / actin gene copy number ratio",trans="log10")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x = element_blank(),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))
dev.off()
med_act_copies <- median(CT_dt[ Target %like% "Actin", copies])
CT_ratio_dt[ , Norm_copies := ratio * med_act_copies ]
CT_ratio_dt <- merge.data.table(CT_ratio_dt, CT_dt[Target %like% "UV", c("Sample", "copiesperngDNA")], by="Sample", all.x = TRUE)


### Read mapping proportions
#Read data
readnum_dt <- fread(paste0(datapath_raw,"readnumbers.csv"))
setnames(readnum_dt, c("Sample","Raw","Bacterial_db","A_mellifera"))
readnum_dt[, Leftover := Raw-Bacterial_db-A_mellifera]
readnum_dt[, prop_bac := Bacterial_db/Raw]
readnum_dt[, prop_amel := A_mellifera/Raw]
reads_abs_stat <- wilcox.test(readnum_dt[Sample %like% "N", prop_bac], 
                              readnum_dt[Sample %like% "F", prop_bac], paired = TRUE)
readnum_mlt <- melt.data.table(readnum_dt, id.vars = "Sample", value.name = "Reads", variable.name = "Mapping")[!(Mapping %like% "prop" | Mapping %like% "Raw"),]

pdf(paste0(figpath, "readnumbers.pdf"), onefile=TRUE, width=13)
ggplot(readnum_mlt[ Mapping!="Raw",],
       aes(x=factor(Sample,
                    levels=sample_order),
           y=Reads,
           fill=factor(Mapping,
                       levels=c("Leftover","A_mellifera","Bacterial_db"),
                       labels=c("Not mapped","A. mellifera genome","Bacterial database")))) +
  geom_bar(stat="identity")+
  xlab("Sample")+
  ylab("Number of reads")+
  scale_fill_discrete(name="Mapping")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5))
dev.off()

read_prop_mlt <- melt.data.table(readnum_dt, id.vars = "Sample", value.name = "Proportion", variable.name = "Mapping")[Mapping %like% "prop",]
pdf(paste0(figpath, "readproportions.pdf"), onefile=TRUE, width=13)
ggplot(read_prop_mlt,
       aes(x=factor(Sample,
                    levels=sample_order),
           y=Proportion,
           fill=factor(Mapping,
                       levels=c("left_prop","prop_amel","prop_bac"),
                       labels=c("Not mapped","A. mellifera genome","Bacterial database")))) +
  geom_bar(stat="identity") +
  xlab("Sample")+
  ylab("Read proportion")+
  scale_fill_discrete(name="Mapping")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5))
dev.off()

read_hostbac_ratio <- copy(readnum_dt)
read_hostbac_ratio[, readratio:=Bacterial_db/A_mellifera]
tmp_ct_reads <- merge.data.table(CT_ratio_dt, read_hostbac_ratio[, c("Sample", "readratio")], by = "Sample")
tmp_ct_reads[ Sample %like% "N", Host := "Nurses" ]
tmp_ct_reads[ Sample %like% "F", Host := "Foragers" ]
ggplot(tmp_ct_reads, aes(x=ratio, y=readratio, fill=factor(Host, levels=treatments)))+
  geom_point(aes(fill = factor(Host, levels = treatments)), 
             shape = 22, size = 6, position = position_dodge2(width = 0.6))+
  labs(subtitle = paste0("Pearson correlation coefficient: \n r = ", 
                         cor(tmp_ct_reads[, ratio], tmp_ct_reads[, readratio], method = "pearson")),
       fill="Host")+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  xlab("16S/actin gene copy number ratio")+
  ylab("bacterial/host read number")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))


### Phylotype relative abundances
#Read data
phylo_all_dt <- data.table()
for (i in phylo_list) {
  tmp_dt <- fread(paste0(datapath, i, "_corecov_coord.txt"), header = TRUE)
  tmp_all_dt <- as.data.table(aggregate(data=tmp_dt, Cov_ter ~ Sample + Cluster, sum))
  tmp_all_dt[, Phylo := i]
  phylo_all_dt <- rbind(phylo_all_dt, tmp_all_dt)
}

#Compute the proportions
cov_sum_dt <- phylo_all_dt[, sum(Cov_ter), by=c("Phylo", "Sample")]
phylo_all_dt <- merge.data.table(phylo_all_dt, cov_sum_dt, by=c("Phylo", "Sample"), all.x = TRUE)
cov_sum_sample_dt <- phylo_all_dt[, sum(Cov_ter), by=Sample]
phylo_all_dt <- merge.data.table(phylo_all_dt, cov_sum_sample_dt, by="Sample", all.x = TRUE)
setnames(phylo_all_dt, c("Sample", "Phylo", "SDP", "Cov_ter", "Phylo_sum", "Sample_sum"))
phylo_all_dt[ , Prop := Cov_ter / Sample_sum ]
phylo_all_dt[ , Prop_SDP := Cov_ter / Phylo_sum]
phylo_all_dt[is.na(Prop_SDP), Prop_SDP:=0]
phylo_all_dt[ , Norm_abun_SDP := setkey(phylo_all_dt, Sample)[ CT_ratio_dt, Prop * copiesperngDNA]]
pdf(paste0(figpath, "phylorelabundstackbarplot.pdf"), onefile=TRUE, width=12)
ggplot(data=phylo_all_dt,
       aes(x=factor(Sample,levels=sample_order),
           y=Prop,
           fill=factor(Phylo,
                       levels=phylo_order_rv,
                       labels=c("Bombella","L. kunkeii","Commensalibacter",
                                "Bartonella","Frischella","Apibacter",
                                "Snodgrasella","Gilliamella","Firm4","Bifidobacterium", "Firm5")))) + 
  geom_bar(stat="identity",colour="white") +
  scale_fill_manual(values=phylo_col_rv, name="Phylotype") +
  xlab("Sample")+
  ylab("Proportion")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5))
dev.off()

#Normalized absolute abundance
phylo_absabun_dt <- as.data.table(aggregate(data=phylo_all_dt, Norm_abun_SDP ~ Sample + Phylo, sum))
setnames(phylo_absabun_dt, c("Sample", "Phylo", "Norm_abun"))
phylo_absabun_dt[ Sample %like% "F[0-9][0-9]", Host:="Foragers"]
phylo_absabun_dt[ Sample %like% "N[0-9][0-9]", Host:="Nurses"]
stats_phylo_absabun <- data.table(phylo_list, 0)
setnames(stats_phylo_absabun, c("Phylo", "p"))
for (i in phylo_list) {
  temp_test <- wilcox.test(phylo_absabun_dt[(Phylo==i & Sample %like% "N[0-9][0-9]"), Norm_abun],
                           phylo_absabun_dt[(Phylo==i & Sample %like% "F[0-9][0-9]"), Norm_abun],
                           paired = TRUE)
  stats_phylo_absabun[ Phylo==i, p:=temp_test$p.value]
}
stats_phylo_absabun[, q:=qvalue(stats_phylo_absabun$p)$qvalues]
phylo_absabun_bplots <- list()
for (i in phylo_list) {
  phylo_absabun_bplots[[i]] <- ggplot(phylo_absabun_dt[ Phylo==i, ],
                                      aes(x=factor(Host,levels=treatments),
                                          y=Norm_abun,
                                          fill=factor(Host,levels=treatments)))+
    geom_boxplot(width=0.5)+
    geom_point(aes(fill = Host),
               shape=22,
               size=4,
               position=position_dodge2(width=0.6)) +
    scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                            "Foragers"="#56B4E9"))+
    ggtitle(i)+
    labs(fill="Host",subtitle=paste0("qvalue = ",stats_phylo_absabun[Phylo==i, q]))+
    scale_y_continuous("Abundance \n (normalized estimated \n 16S gene copies)",
                       trans="sqrt")+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          text=element_text(size=22),
          axis.title.x=element_blank(),
          plot.title=element_text(hjust=0.5), 
          plot.subtitle=element_text(hjust=0.5))
}
pdf(paste0(figpath, "phylorabsabundboxplot.pdf"), onefile=TRUE, width=12)
marrangeGrob(phylo_absabun_bplots, nrow = 1, ncol = 1)
dev.off()

### SDP relative abundances and Nurses vs. Foragers normalized abundances boxplots
stats_SDP_absabun <- data.table(unique(phylo_all_dt[, SDP]), 0)
setnames(stats_SDP_absabun, c("SDP", "p"))
for (i in unique(phylo_all_dt[, SDP])) {
  temp_test <- wilcox.test(phylo_all_dt[(SDP==i & Sample %like% "N[0-9][0-9]"), Norm_abun_SDP],
                           phylo_all_dt[(SDP==i & Sample %like% "F[0-9][0-9]"), Norm_abun_SDP],
                           paired = TRUE)
  stats_SDP_absabun[ SDP==i, p:=temp_test$p.value]
}
stats_SDP_absabun[, q:=qvalue(stats_SDP_absabun$p)$qvalues]
phylo_all_dt[ Sample %like% "F[0-9][0-9]", Host:="Foragers"]
phylo_all_dt[ Sample %like% "N[0-9][0-9]", Host:="Nurses"]
sdp_relabun_barplot <- list()
sdp_absabun_boxplot <- list()
for (i in phylo_list) {
  #Relative abundances stuff
  sdp_relabun_barplot[[i]] <- ggplot(data=phylo_all_dt[Phylo==i,],
                                     aes(x=factor(Sample,levels=sample_order),
                                         y=Prop_SDP,
                                         fill=SDP)) + 
    geom_bar(stat="identity", colour="white") +
    scale_fill_manual(values=phylo_col_rv, name = "SDP") + 
    ggtitle(paste0(i," SDP relative abundances")) +
    xlab("Sample")+
    ylab("Proportion")+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          text=element_text(size=22),
          plot.title=element_text(hjust=0.5))
}
pdf(paste0(figpath, "sdprelabunstackbarplot.pdf"), onefile=TRUE, width = 12)
marrangeGrob(sdp_relabun_barplot, nrow = 1, ncol = 1)
dev.off()
stats_SDP_absabun <- merge.data.table(stats_SDP_absabun,
                                      phylo_all_dt[Sample=="F01", c("SDP","Phylo")], 
                                      by="SDP", all.x=TRUE )
pdf(paste0(figpath, "sdpabsabunboxplot.pdf"), onefile=TRUE, width = 12)
ggplot(phylo_all_dt[!(SDP %in% c("api_apis_dorsa", "api_bombus", "bifido_1_cerana", 
                                 "bifido_bombus","bom_1", "bom_apis_melli","bom_bombus",
                                 "com_drosophila", "com_monarch", "firm5_bombus",
                                 "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre",
                                 "gilli_apis_dorsa", "gilli_bombus","snod_2", "snod_bombus")),], 
       aes(x=SDP,
           y=Norm_abun_SDP))+
  geom_boxplot(aes(fill=factor(Host, levels=treatments)))+
  geom_point(data=stats_SDP_absabun[q<0.05 & !(SDP %in% c("api_apis_dorsa", "api_bombus", "bifido_1_cerana", 
                                                          "bifido_bombus","bom_1", "bom_apis_melli","bom_bombus",
                                                          "com_drosophila", "com_monarch", "firm5_bombus",
                                                          "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre",
                                                          "gilli_apis_dorsa", "gilli_bombus","snod_2", "snod_bombus")), ], 
             aes(x=SDP, y=135000), shape="*", size=8, show.legend = FALSE)+
  xlab("SDP")+
  scale_fill_manual(name="Host", values=c("Nurses"="#D55E00",
                                          "Foragers"="#56B4E9"))+
  facet_grid(. ~ Phylo, scales="free", space="free")+
  scale_y_continuous(name="Abundance \n (normalized estimated \n 16S gene copies)", trans="sqrt")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
dev.off()

phylo_all_dt[, Hive:=gsub("[FN]","", Sample)]
sdp_fc_dt <- dcast.data.table(phylo_all_dt, Hive + SDP + Phylo ~ Host, value.var = "Norm_abun_SDP", )
sdp_fc_dt[Nurses>0, logtwoFC:=log2(Foragers/Nurses)]

sdp_fc_dt_bis <- copy(sdp_fc_dt)
sdp_fc_dt_bis <- sdp_fc_dt_bis[ !(SDP %in% c("api_apis_dorsa", "api_bombus", "bifido_1_cerana", 
                            "bifido_bombus", "bom_apis_melli","bom_bombus",
                            "com_drosophila", "com_monarch", "firm5_bombus",
                            "gilli_4", "gilli_5", "gilli_6", "gilli_apis_andre",
                            "gilli_apis_dorsa", "gilli_bombus", "snod_bombus")),]

pdf(paste0(figpath, "fc_sdp_plot.pdf"), onefile=TRUE, width=15)
ggplot(sdp_fc_dt_bis[!(SDP %in% c("api_1", "bom_1", "firm5_7", "snod_2") ),], 
       aes(x = SDP, 
           y = logtwoFC))+
  geom_point(aes(fill=factor(Phylo,
                             levels=phylo_order_rv,
                             labels = c("Bombella", "L. kunkeii", "Commensalibacter",
                                        "Bartonella", "Frischella", "Apibacter",
                                        "Snodgrasella", "Gilliamella", "Firm4", "Bifidobacterium", "Firm5"))), 
             shape = 22, size = 3, position = position_dodge2(width = 0.2)) +
  geom_hline(yintercept = 0)+
  geom_crossbar(data=sdp_fc_dt_bis[!is.na(logtwoFC) & !(SDP %in% c("api_1", "bom_1", "firm5_7", "snod_2")), median(logtwoFC), by=SDP],
                aes(x=SDP,ymin=V1, ymax=V1,y=V1,group=SDP), width = 0.5)+
  geom_point(data=stats_SDP_absabun[q<0.05 & !(SDP %in% c("api_1", "bom_1", "firm5_7", "snod_2")), ], 
             aes(x=SDP, y=3), shape="*", size=8, show.legend = FALSE)+
  scale_fill_manual(values=phylo_col_rv[c(2,3,4,5,7,8,9,10,11)], name = "Phylotype") + 
  xlab("SDP")+
  scale_y_continuous("log2 fold change to nurse sample")+
  theme(axis.text.x=element_text(angle=45,hjust=1), 
        text=element_text(size=22),
        plot.title=element_text(hjust=0.5), 
        plot.subtitle=element_text(hjust=0.5))
dev.off()

### NMDS plots for the relative and absolute compositions of samples, based on phylotypes and based on SDPs, using Bray-Curtis dissimilarity indices.

nmds_plot <- function(data_dt, plot_title){
  set.seed(123)
  nmds = metaMDS(data_dt, distance = "bray")
  nmds_scores = as.data.table(scores(nmds))
  nmds_scores[, Sample := phylo_dt_umlt_prop$Sample]
  nmds_scores[Sample %like% "N[0-9][0-9]", Host := "Nurses" ]
  nmds_scores[Sample %like% "F[0-9][0-9]", Host := "Foragers" ]
  nmds_plot <- ggplot(nmds_scores, aes(x=NMDS1, y=NMDS2, 
                                     color=factor(Host,levels=treatments),
                                     shape=factor(Host,levels=treatments)))+
    geom_point(stat="identity",size=6)+
    labs(x="NMDS1",
         y="NMDS2",
         title=plot_title,
         subtitle=paste0("Stress = ",nmds$stress),
         color="Host",
         shape="Host")+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=22),
          plot.subtitle = element_text(hjust = 0.5))
  return(nmds_plot)
}

#Using phylotype proportions
phylo_dt_umlt_prop <- dcast.data.table(phylo_all_dt, Sample ~ Phylo, value.var = "Prop", fun.aggregate = sum)
phylo_prop_mat <- as.matrix(phylo_dt_umlt_prop[, -1])
phylo_prop_nmdsplot <- nmds_plot(phylo_prop_mat, "Phylotype proportions-based NMDS plot")

#Using phylotype absolute abundances
phylo_dt_umlt_absabun <- dcast.data.table(phylo_absabun_dt, Sample ~ Phylo, value.var = "Norm_abun")
phylo_absabun_mat <- as.matrix(phylo_dt_umlt_absabun[, -1])
phylo_absabun_nmdsplot <- nmds_plot(phylo_absabun_mat, "Phylotype normalized absolute abundances-based NMDS plot")

#Using SDP proportions
sdp_dt_umlt_prop <- dcast.data.table(phylo_all_dt, Sample ~ SDP, value.var = "Prop")
sdp_prop_mat <- as.matrix(sdp_dt_umlt_prop[,-1])
sdp_prop_nmdsplot <- nmds_plot(sdp_prop_mat, "SDP proportions-based NMDS plot")

#Using SDP absolute abundances
sdp_dt_umlt_absabun <- dcast.data.table(phylo_all_dt, Sample ~ SDP, value.var = "Norm_abun_SDP")
sdp_absabun_mat <- as.matrix(sdp_dt_umlt_absabun[,-1])
sdp_absabun_nmdsplot <- nmds_plot(sdp_absabun_mat, "SDP normalized absolute abundances-based NMDS plot")

pdf(paste0(figpath, "phyloandsdp_NDMS.pdf"), onefile=TRUE, width = 12)
marrangeGrob(list(phylo_prop_nmdsplot, phylo_absabun_nmdsplot, sdp_prop_nmdsplot, sdp_absabun_nmdsplot), nrow = 1, ncol = 1)
dev.off()


pdf(paste0(figpath, "nmds_rel_phylo.pdf"), onefile=TRUE, width=10)
phylo_prop_nmdsplot
dev.off()
pdf(paste0(figpath, "nmds_rel_sdp.pdf"), onefile=TRUE, width=10)
sdp_prop_nmdsplot
dev.off()

