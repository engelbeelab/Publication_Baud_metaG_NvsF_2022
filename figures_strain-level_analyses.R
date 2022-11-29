library(ggplot2)
library(ape)
library(gridExtra)
library(RColorBrewer)
library(reshape)
library(ggrepel)
library(qvalue)
library(data.table)
library(scales)
library(PERMANOVA)

#Filepaths
qPCR_datapath <- "./"
datapath <- "./strain-level_analysis/"
figpath <- "./strain-level_analysis/"
metadatapath <- "./db_metadata_GB.txt"

### Sample lists
meta_dt <- fread(metadatapath, header = TRUE)
Nlist <- c("N01","N02","N03","N04","N06","N07","N08","N09","N10","N11","N12","N13","N14","N15","N16")
Flist <- c("F01","F02","F03","F04","F06","F07","F08","F09","F10","F11","F12","F13","F14","F15","F16")
samplelist <- c(Nlist, Flist)
snod_1_Nlike_strains <- c("Ga0326500","CPT77","BGH94","BGH95","BGH96","BGH97","BGI00", 
                          "BGI01","BGI08","BGI09","BGI12","BGI13","BHC42","BHC50", 
                          "BHC52","Ga0227306","SALWKB2","H3V01","H3V10","H3V11",
                          "H3U73","H3U82","H3T75","Ga0227304")
snod_1_Flike_strains <- c("BHC47","Ga0227305","H3T79")
sample_order <- c(Nlist, Flist)
phylo_order <- c("firm5","bifido","firm4","gilli","snod","api","fper","bapis","com","lkun", "bom")
phylo_order_rv <- rev(phylo_order)
phylo_list <- c("api","bapis","bifido","bom","com","firm4","firm5","fper","gilli","lkun","snod")
sdp_list <- c("api_1","bapis","bifido_1.1","bifido_1.2",
              "bifido_1.3","bifido_1.4","bifido_1.5","bifido_2","bom_1","com_1",
              "firm4_1","firm4_2","firm5_1","firm5_2","firm5_3","firm5_4","firm5_7",
              "fper_1","gilli_1","gilli_2","gilli_3","gilli_4","gilli_5","gilli_6",
              "lkun","snod_1","snod_2")
sdp_list_restricted <- c("bapis","bifido_1.1","bifido_1.2","bifido_1.3","bifido_1.4",
                         "bifido_1.5","bifido_2","firm4_1","firm4_2","firm5_1",
                         "firm5_2","firm5_3","firm5_4","fper_1","gilli_1",
                         "gilli_2","snod_1")
treatments <- c("Nurses","Foragers")
treat_colors <- c("#D55E00","#56B4E9","#F0E442")

### qPCR data import
#qPCR standard curve values from Kesnerova et al. (2017)
actin_intercept = 35.0942119870867
actin_slope = -3.2699442250388
UV_intercept = 36.5821936471122
UV_slope = -3.35085896083287

filepath <- paste0(qPCR_datapath, "20190813_AllHives.csv")
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
CT_ratio_dt <- merge.data.table(CT_dt[ Target %like% "Actin",],
                                CT_dt[ Target %like% "UV",],
                                by=c("Sample"))
CT_ratio_dt[, ratio := copies.y / copies.x]
med_act_copies <- median(CT_dt[ Target %like% "Actin", copies])
CT_ratio_dt[ , Norm_copies := ratio * med_act_copies ]
CT_ratio_dt <- merge.data.table(CT_ratio_dt, CT_dt[Target %like% "UV", c("Sample", "copiesperngDNA")], by="Sample", all.x = TRUE)

### Computing total fraction of SNVs per SDP, per host

var_dt <- fread(paste0(datapath, "tot_var_all.txt"),h=F)
setnames(var_dt, c("Phylo","SDP","Host","Nb_snv","Fraction_var","Nb_samples","Cum_ter_cov"))

pdf(paste0(figpath, "snvs_tot_var_all.pdf"), onefile=TRUE, width=12)
ggplot(var_dt[Host!="g",], aes(x=Nb_samples,
                               y=Fraction_var, 
                               color=factor(Host,
                                            levels=c("N","F"),
                                            labels=treatments)))+
  geom_point(aes(size=Cum_ter_cov))+
  scale_size(name="Cumulative terminus coverage")+
  xlab("# samples") +
  ylab("Total polymorphic sites per SDP (%)")+
  scale_color_manual(values=treat_colors,name="Host")+
  theme(text=element_text(size=16),
        panel.background=element_rect(fill="light grey"),
        panel.grid.major=element_line(color="white"))
dev.off()

### Computing fractions of SNVs per SDP, per sample

frac_snvs_dt <- fread(paste0(datapath, "all_sample_var_host.txt"),h=F)
setnames(frac_snvs_dt, c("Phylo", "SDP","Host","Colony","Sample", "Nb_snps","Fraction_var"))
aggr <- as.data.table(aggregate(Fraction_var ~ SDP + Host + Phylo, frac_snvs_dt, median))

stats_frac_var <- data.table(sdp_list_restricted, 0)
setnames(stats_frac_var, c("SDP", "p"))
for (i in sdp_list_restricted) {
  temp_test <- wilcox.test(frac_snvs_dt[(SDP==i & Sample %like% "N[0-9][0-9]"), Fraction_var],
                           frac_snvs_dt[(SDP==i & Sample %like% "F[0-9][0-9]"), Fraction_var],
                           paired = FALSE)
  stats_frac_var[ SDP==i, p:=temp_test$p.value]
}
stats_frac_var[, q:=qvalue(stats_frac_var$p)$qvalues]
stats_frac_var <- merge.data.table(stats_frac_var,
                                      frac_snvs_dt[, .N, by=c("SDP","Phylo")][, c("SDP", "Phylo")], 
                                      by="SDP", all.x=TRUE )

frac_snvs_dt <- frac_snvs_dt[!(Host=="g"),]

pdf(paste0(figpath, "var_per_sample.pdf"), onefile=TRUE, width=12)
ggplot(data=frac_snvs_dt[SDP %in% sdp_list_restricted,], 
       aes(x=SDP,
           y=Fraction_var)) +
  geom_jitter(position=position_jitterdodge(), alpha=0.5,
              aes(color=factor(Host,levels=c("N","F"),labels=treatments))) +
  geom_point(data = aggr[!(Fraction_var)==0 & SDP %in% sdp_list_restricted,], 
             aes(y=Fraction_var, 
                 x=SDP,
                 fill=factor(Host,levels=c("N","F"),labels=treatments),
                 color=factor(Host,levels=c("N","F"),labels=treatments)),
             position=position_jitterdodge(),
             colour="black",
             show.legend = FALSE,
             shape=3,
             alpha=1)+
  geom_point(data=stats_frac_var[q<0.05, ],
             aes(x=SDP, y=8.5), shape="*", size=8, show.legend = FALSE)+
  theme(axis.text.x=element_text(angle=90),
        text=element_text(size = 16), 
        panel.background=element_rect(fill="light grey"), 
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white")) +
  xlab("SDP")+
  scale_y_continuous(name="Polymorphic sites per sample (%)", limits=c(0,8.5), n.breaks=9)+
  # facet_grid(Phylo ~ ., scales="free")+
  facet_grid(. ~ Phylo, scales="free", space="free")+
  scale_color_manual(values=treat_colors, name="Host")#+
  #coord_flip()
dev.off()

### SNV correlations with ct ratios
stats_corr_ct_snvs <- data.table(sdp_list_restricted, 0,0,0,0)
setnames(stats_corr_ct_snvs, c("SDP", "rN","rF","pN","pF"))
tmp_dt <- copy(frac_snvs_dt)
tmp_m_dt <- merge.data.table(tmp_dt, CT_dt[Target %like% "UV", c("Sample", "copiesperngDNA")], by="Sample")
tmp_phylo_all_dt <- copy(phylo_all_dt)
tmp_phylo_all_dt[SDP=="fper", SDP:="fper_1"]
tmp_m_dt <- merge.data.table(tmp_m_dt, tmp_phylo_all_dt[, c("Sample","SDP", "Cov_ter", "Prop")], by=c("Sample","SDP"), all.x=TRUE)
tmp_m_dt[, ratio := Cov_ter/(copiesperngDNA*Prop)]
tmp_m_dt[Host=="N", Host:="Nurses"]
tmp_m_dt[Host=="F", Host:="Foragers"]
for(i in sdp_list_restricted){
  rn <- format(round(cor(tmp_m_dt[SDP==i & Host=="Nurses", ratio],
                         tmp_m_dt[SDP==i & Host=="Nurses", Fraction_var], method = "pearson"),3),nsmall=3)
  stats_corr_ct_snvs[SDP==i, rN:=rn]
  rf <- format(round(cor(tmp_m_dt[SDP==i & Host=="Foragers", ratio],
                         tmp_m_dt[SDP==i & Host=="Foragers", Fraction_var], method = "pearson"),3),nsmall=3)
  stats_corr_ct_snvs[SDP==i, rF:=rf]
  pn <- try(format(round(cor.test(tmp_m_dt[SDP==i & Host=="Nurses", ratio],
                                  tmp_m_dt[SDP==i & Host=="Nurses", Fraction_var], method = "pearson")$p.value,3),
                   nsmall=3),silent=TRUE)
  stats_corr_ct_snvs[SDP==i, pN:=pn]
  pf <- try(format(round(cor.test(tmp_m_dt[SDP==i & Host=="Foragers", ratio],
                                  tmp_m_dt[SDP==i & Host=="Foragers", Fraction_var], method = "pearson")$p.value,3),
                   nsmall=3),silent=TRUE)
  stats_corr_ct_snvs[SDP==i, pF:=pf]
}
stats_corr_ct_snvs <- merge.data.table(stats_corr_ct_snvs, tmp_m_dt[, c("SDP","Phylo")], by="SDP", all.x=TRUE)
pdf(paste0(figpath, "snvs_ctratio_correlation2.pdf"), onefile=TRUE, width=12, height=30)
ggplot(tmp_m_dt[SDP %in% sdp_list_restricted,], aes(x=as.numeric(ratio), y=Fraction_var))+
  geom_point(aes(fill=factor(Host,levels=treatments),
                 shape=factor(Host,levels=treatments)), 
             size=2)+
  scale_shape_manual(name="Host",values=c(21,24))+
  scale_fill_manual(name="Host", values=treat_colors)+
  geom_text(data=stats_corr_ct_snvs, 
            aes(x = Inf, y = -Inf, label=paste0("r(Nurses)=",rN,"\n",
                                                "p(Nurses)=",pN,"\n",
                                                "r(Foragers)=",rF,"\n",
                                                "p(Foragers)=",pF)), 
            show.legend = FALSE,
            vjust=-0.5,
            hjust=1.1,
            size=4)+
  xlab("Terminus coverage (reads/bp) \n per bacterial load (SDP proportion * 16S copy/ng of DNA)")+
  ylab("Polymorphic sites per sample (%)")+
  facet_grid(SDP ~ ., scales="free_y", space="fixed")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16),
        panel.background = element_rect(fill = "light grey"), 
        panel.grid.major = element_line(color = "white"))
dev.off()

### SNVs cumulative curves
cum_curves_dt <- data.table()
for (i in sdp_list_restricted) {
  tmp <-  fread(paste0(datapath,i, "_cum_curve.txt"), h=F)
  cum_curves_dt <- rbind(cum_curves_dt, tmp)
}
setnames(cum_curves_dt, c("Host","curve","Nsamples","percvar","SDP"))
pdf(paste0(figpath, "snvs_cum_curves.pdf"), onefile=TRUE, width=12)
ggplot(data=cum_curves_dt[(Host=="F" | Host=="N"),], 
       aes(x=Nsamples, 
           y=percvar,
           colour=factor(Host,
                         levels=c("N", "F"),
                         labels=treatments))) +
  geom_jitter(position=position_dodge(width=0.7)) + 
  geom_smooth(se=FALSE)+
  xlab("# Bee samples") + 
  ylab("Polymorphic sites (%)") +
  scale_x_continuous(breaks=c(0,15))+
  scale_color_manual(values = treat_colors, name="Host")+
  facet_wrap(. ~ SDP, ncol=6,nrow=3)+
  theme(panel.background=element_rect(fill="light grey"),
        panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"), 
        panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"), 
        legend.title=element_blank())
dev.off()


### PCoAs

#pcoa plot function
pcoa_plot<- function(dt, graph_title, locations_enabled=TRUE) {
  dt_c <- copy(dt)
  db_cols <- dt_c[!(V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]"), V1]
  dt_c <- dt_c[V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]", ]
  dt_c[, (db_cols):=NULL]

  if (nrow(dt_c[V1 %like% "F[0-9][0-9]",])==0 | nrow(dt_c[V1 %like% "N[0-9][0-9]",]) == 0) {
    return()
  }
  matrix <- as.matrix(dt_c, rownames = "V1")
  dist <- as.dist(matrix)
  res_pcoa <- pcoa(dist)
  ev1 <- res_pcoa$vectors[,1]
  ev2 <- res_pcoa$vectors[,2]
  df_new <- data.frame(cbind(ev1,ev2))    
  sample_groups <- substr(rownames(df_new),1,1)
  df_new <- cbind(df_new,sample_groups)
  
  Di <- list()
  Di[["D"]] <- matrix
  fac<-as.factor(substr(dt_c$V1,1,1))
  perm <- PERMANOVA(Di, fac)
  Fstat <- eval(parse(text=perm$Initial$Global)[5])
  
  locations <- rownames(df_new)
  hive <- gsub("[N F]", "", locations)
  df_new$hive <- hive
  locations <- gsub("^[A-E G-M O-Z].*", "dbstrains", locations)
  locations <- gsub("[N F]0[1-3]", "UNIL", locations)
  locations <- gsub("[N F]0[4-7]", "Liebefeld", locations)
  locations <- gsub("[N F]0[8-9]", "Yens", locations)
  locations <- gsub("[N F]10", "Yens", locations)
  locations <- gsub("[N F]1[1-3]", "Cugy", locations)
  locations <- gsub("[N F]1[4-6]", "Vesancy", locations)
  df_new$locations <- locations
  
  perc_axis <- round(((res_pcoa$values$Relative_eig[c(1,2)])*100), digits=1)
  axis_x_title <- paste0("PCo1 (",perc_axis[1],"%)")
  axis_y_title <- paste0("PCo2 (",perc_axis[2],"%)")
  
  if (!locations_enabled) {
    p <- ggplot(df_new,aes(x=ev1,
                           y=ev2,
                           color=factor(sample_groups,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers")), 
                           shape=factor(sample_groups,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers"))))+
      geom_point(stat="identity",size=2)+
      labs(x=as.character(axis_x_title),
           y=as.character(axis_y_title), 
           title=graph_title,
           subtitle = paste0("Permanova: \n Pseudo-F = ",Fstat))+
      scale_color_manual(values = treat_colors, name = "Host")+
      scale_shape_discrete(name="Host")+
      theme(legend.title =element_blank(),
            panel.background=element_rect(fill="light grey"),
            panel.grid.major=element_line(size=0.5,linetype='solid',colour="white"), 
            panel.grid.minor=element_line(size=0.25,linetype='solid',colour="white"),
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8),
            plot.title=element_text(size=10,face="plain"),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8))
  } else {
    p <- ggplot(df_new,aes(x=ev1,
                           y=ev2,
                           colour=factor(locations), 
                           alpha=factor(sample_groups,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers")),
                           shape=factor(sample_groups,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers")),
                           size=factor(sample_groups,
                                        levels=c("N","F"),
                                        labels=c("Nurses","Foragers"))))+
      geom_point(stat="identity")+
      geom_line(aes(group=factor(hive)), alpha=0.5, size=0.5)+
      scale_size_manual(values=c("Foragers"=2,"Nurses"=3), guide="none")+
      scale_alpha_manual(values = c("Foragers"=1, "Nurses"=0.5), guide='none')+
      scale_color_manual(values=c("UNIL"="#CC79A7",
                                  "Cugy"="#0072B2",
                                  "Yens"="#009E73",
                                  "Liebefeld"="#E69F00",
                                  "Vesancy"="#999999"),
                         name = "Host")+
      labs(x=as.character(axis_x_title),
           y=as.character(axis_y_title), 
           title=graph_title)+
      theme(legend.title =element_blank(),
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8),
            plot.title=element_text(size=10,face="plain"),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8))
  }
  
  return(p)
}
get_permanova<- function(dt) {
  dt_c <- copy(dt)
  db_cols <- dt_c[!(V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]"), V1]
  dt_c <- dt_c[V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]", ]
  dt_c[, (db_cols):=NULL]
  
  if (nrow(dt_c[V1 %like% "F[0-9][0-9]",])==0 | nrow(dt_c[V1 %like% "N[0-9][0-9]",]) == 0) {
    return()
  }
  matrix <- as.matrix(dt_c, rownames = "V1")
  dist <- as.dist(matrix)
  res_pcoa <- pcoa(dist)
  ev1 <- res_pcoa$vectors[,1]
  ev2 <- res_pcoa$vectors[,2]
  df_new <- data.frame(cbind(ev1,ev2))    
  sample_groups <- substr(rownames(df_new),1,1)
  df_new <- cbind(df_new,sample_groups)
  
  Di <- list()
  Di[["D"]] <- matrix
  fac<-as.factor(substr(dt_c$V1,1,1))
  perm <- PERMANOVA(Di, fac)
  return(perm)
}
get_permanova_loc <- function(dt) {
  dt_c <- copy(dt)
  db_cols <- dt_c[!(V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]"), V1]
  dt_c <- dt_c[V1 %like% "F[0-9][0-9]" | V1 %like% "N[0-9][0-9]", ]
  dt_c[, (db_cols):=NULL]
  
  if (nrow(dt_c[V1 %like% "F[0-9][0-9]",])==0 | nrow(dt_c[V1 %like% "N[0-9][0-9]",]) == 0) {
    return()
  }
  matrix <- as.matrix(dt_c, rownames = "V1")
  dist <- as.dist(matrix)
  res_pcoa <- pcoa(dist)
  ev1 <- res_pcoa$vectors[,1]
  ev2 <- res_pcoa$vectors[,2]
  df_new <- data.frame(cbind(ev1,ev2))    
  sample_groups <- rownames(df_new)
  sample_groups <- gsub("[N F]0[1-3]", "UNIL", sample_groups)
  sample_groups <- gsub("[N F]0[4-7]", "Liebefeld", sample_groups)
  sample_groups <- gsub("[N F]0[8-9]", "Yens", sample_groups)
  sample_groups <- gsub("[N F]10", "Yens", sample_groups)
  sample_groups <- gsub("[N F]1[1-3]", "Cugy", sample_groups)
  sample_groups <- gsub("[N F]1[4-6]", "Vesancy", sample_groups)
  df_new <- cbind(df_new,sample_groups)
  
  Di <- list()
  Di[["D"]] <- matrix
  fac<-as.factor(substr(df_new$sample_groups,1,1))
  perm <- PERMANOVA(Di, fac)
  return(perm)
}

#Reading the data
snvs_sdp_dt <- list()
snvs_sdp_plots <- list()
snvs_sdp_plots_loc <- list()
permanovastats <- list()
permanova_vals <- data.table()
permanova_locstats <- list()

#Making the plots
for (i in sdp_list_restricted) {
  snvs_sdp_dt[[i]] <- fread(paste0(datapath, i, "_dist_matrix_KE.txt"),h=T)
  snvs_sdp_dt[[i]] <- snvs_sdp_dt[[i]][, c("V1", snvs_sdp_dt[[i]][, V1]), with=FALSE]
  snvs_sdp_plots[[i]] <- pcoa_plot(snvs_sdp_dt[[i]], i, locations_enabled = FALSE)
  snvs_sdp_plots_loc[[i]] <- pcoa_plot(snvs_sdp_dt[[i]], i, locations_enabled = TRUE)
  permanovastats[[i]] <- get_permanova(snvs_sdp_dt[[i]])
  permanova_locstats[[i]] <- get_permanova_loc(snvs_sdp_dt[[i]])
  permanova_vals <- rbind(permanova_vals, data.table(i,eval(parse(text=permanovastats[[i]]$Initial$Global)[6]),eval(parse(text=permanovastats[[i]]$Initial$Global)[5]),
                                                     eval(parse(text=permanova_locstats[[i]]$Initial$Global)[6]),eval(parse(text=permanova_locstats[[i]]$Initial$Global)[5])))
}
setnames(permanova_vals, c("SDP","pvalue","Fstat","locpval","locFstat"))
permanova_vals[, qval:=qvalue(as.numeric(permanova_vals[, pvalue]), pi0 = 1)$q]
permanova_vals[, locqval:=qvalue(as.numeric(permanova_vals[, locpval]), pi0 = 1)$q]
pdf(paste0(figpath, "snvs_pcoa_nolabel.pdf"), onefile=TRUE,width=15,height=15)
marrangeGrob(snvs_sdp_plots, nrow = 5, ncol = 4, layout_matrix = matrix(1:20,5,4,TRUE))
dev.off()
pdf(paste0(figpath, "snvs_pcoa_loc.pdf"), onefile=TRUE,width=15,height=15)
marrangeGrob(snvs_sdp_plots_loc, nrow = 5, ncol = 4, layout_matrix = matrix(1:20,5,4,TRUE))
dev.off()

### Boxplots 
snvs_counts_list <- list()
snvs_counts_all <- data.table()
for (i in sdp_list_restricted) {
  snvs_counts_list[[i]] <- fread(paste0(datapath, i, "_filt_snvs_diff_count.txt"), header=FALSE)
  tmp_dt <- copy(snvs_counts_list[[i]])
  tmp_dt[V1 %in% Nlist & V2 %in% Nlist, comparison:="NN"]
  tmp_dt[V1 %in% Flist & V2 %in% Flist, comparison:="FF"]
  tmp_dt[!(V1 %in% samplelist) & !(V2 %in% samplelist), comparison:="gg"]
  tmp_dt[(V1 %in% Nlist & V2 %in% Flist) | (V1 %in% Flist & V2 %in% Nlist), comparison:="NF"]
  tmp_dt[(V1 %in% Nlist & !(V2 %in% samplelist)) | (!(V1 %in% samplelist) & V2 %in% Nlist), comparison:="gN"]
  tmp_dt[(!(V1 %in% samplelist) & V2 %in% Flist) | (V1 %in% Flist & !(V2 %in% samplelist)), comparison:="gF"]
  tmp_dt[, SDP:=i]
  snvs_counts_all <- rbind(snvs_counts_all, tmp_dt)
}
snvs_counts_all[, same_hive:="no"]
snvs_counts_all[comparison=="NF" & substr(V1,2,3)==substr(V2,2,3), same_hive:="yes"]
snvs_counts_all <- merge.data.table(snvs_counts_all,
                                    frac_snvs_dt[, .N, by=c("SDP","Phylo")][, c("SDP", "Phylo")], 
                                    by="SDP", all.x=TRUE )

stats_snvs_counts <- rbind(data.table(sdp_list_restricted, "no", 1),
                           data.table(sdp_list_restricted, "yes", 1))
setnames(stats_snvs_counts, c("SDP", "same_hive", "p"))
for (i in sdp_list_restricted) {
  if(length(snvs_counts_all[(SDP==i & comparison=="NF" & same_hive=="no"), V3])>0 &
     length(snvs_counts_all[(SDP==i & comparison=="NF" & same_hive=="yes"), V3])>0) {
    temp_test <- wilcox.test(snvs_counts_all[(SDP==i & same_hive=="no"), V3],
                             snvs_counts_all[(SDP==i & same_hive=="yes"), V3],
                             paired = FALSE)
    stats_snvs_counts[ SDP==i & same_hive=="yes", p:=temp_test$p.value]
  }
}
stats_snvs_counts[, q:=1]
stats_snvs_counts[same_hive=="yes", q:=qvalue(stats_snvs_counts[same_hive=="yes", p])$qvalues]
stats_snvs_counts <- merge.data.table(stats_snvs_counts,
                                      frac_snvs_dt[, .N, by=c("SDP","Phylo")][, c("SDP", "Phylo")], 
                                      by="SDP", all.x=TRUE )
stats_snvs_counts <- stats_snvs_counts[, min(q), by=c("SDP", "Phylo")] #Significance when both NN v NF and FF v NF comparisons are significantly different.

pdf(paste0(figpath, "snvs_diff_counts_boxplots.pdf"), onefile=TRUE, width = 13)
ggplot(data=snvs_counts_all[comparison=="NF",], aes(x=SDP, y=V3))+
  geom_boxplot(aes(fill=factor(same_hive, 
                               levels=c("no","yes"))))+
  geom_point(data=stats_snvs_counts[V1<0.05,],
             aes(x=SDP, y=60000), shape="*", size=8, show.legend = FALSE)+
  xlab("SDP")+
  scale_fill_manual(values = c("#999999","#F0E442"), name="Same hive of origin?")+
  facet_grid(. ~ Phylo, scales="free", space="free")+
  scale_y_continuous(name="Number of differing SNVs in Forager-Nurse sample pairs", trans="sqrt")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
dev.off()





