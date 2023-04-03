##############
#read in data#
##############

# Liver RNA seq for C57Bl6 mice
load("~/Documents/Work/CBMR/Projects/Seasonal Study/RNA/Liver RNA-seq/Circadian Analysis/Seasonal Liver RNAseq/dataDGE.groups.Rdata")
Pooled.Metadata <- read.csv("~/Documents/Work/CBMR/Projects/Seasonal Study/RNA/Liver RNA-seq/Circadian Analysis/Seasonal Liver RNAseq/Pooled Metadata.csv")

############
#Packages#
###########
library(tidyverse)
library(edgeR)
library(Rmisc)
library(compareRhythms)
library(patchwork)
library(clusterProfiler)
library(reshape2)
library(eulerr)
library(DGEobj.utils)
library(stringr)
library(ggplotify)
library(data.table)

############   
#Functions#
###########

#Panel of clock genes
clock_panel <-  c("Arntl","Nr1d1","Dbp","Npas2","Per1","Per2")

# Main photoperiod colours for figures 
Photocolour <-  c("#56B4E9","#000000", "#FFCC33")

# sine and spline curves for graphing in geom_smooth
sine_curve_fit <-  y~cos(pi*x/12) + sin(pi*x/12)
cosinor_curve <- geom_smooth(
  se = FALSE,
  span = 0.5,
  method = "lm",
  formula = y ~ cos(pi * x / 12) + sin(pi * x / 12))
spline_curve <- geom_smooth(se = FALSE)

#Timer function
# This function converts times in (h) that are negative or over 24 to 24 hour time for subtracting one time from another for phase changes 
timer = function(time){
  if (time > 24)
    timeres = time - 24
  if (time < 0)
    timeres = time + 24
  if (time > 0 & time < 24)
    timeres = time
  return(timeres)
}
timerv= Vectorize(timer)

# Hour calculator function
# This function calculates the smallest time gap between two 24 hour times 
hourCalculator <- function(time1, time2){
  timeResults <- time2 - time1
  if(timeResults>12){
    timeResults <- timeResults - 24
  }
  if(timeResults<(-12)){
    timeResults <- timeResults + 24
  }
  return(timeResults)
}
hourCalculatorv <- Vectorize(hourCalculator)

# This function plots the C57BL6 data in a circadian fashion combining diets with mean ± SE
circadian_plotter_BL6 <- function(gene,sheet,time,columns,smoothing){
  
  d <- filter(summary,{{gene}} %in% sheet)
    graph <- ggplot(d,aes_string(x = time, y = "LogCPM", color = "Photoperiod")) + geom_linerange(aes(ymin=LogCPM, ymax=LogCPM+se)) +geom_point() + smoothing + theme_classic() + facet_wrap( ~ Symbol, scale = "free_y",ncol = columns)+ scale_color_manual(values=Photocolour)+scale_x_continuous(breaks=c(0,6,12,18,24),limits = c(0,24))+theme(strip.text = element_text(face = "italic"))
    return(graph)}

########################
#Voom and limma of Liver RNA seq
########################

# First we use  Voom and limma to look at changes in differential expression between diets and photoperiods regardless of time

# make factors
dataDGE$samples <- dplyr::select(dataDGE$samples,-group)
dataDGE$samples$Photoperiod <- factor(dataDGE$samples$Photoperiod,
                                      levels = c("SL", "EL", "LL"))
dataDGE$samples$Diet <- factor(dataDGE$samples$Diet,
                               levels = c("LFD", "HFD"))
dataDGE$samples$zTimeFact <- as.factor(dataDGE$samples$zTime)

dataDGE$samples$subGroup <- paste(dataDGE$samples$Photoperiod,sep = "_",dataDGE$samples$Diet)
dataDGE$samples$subGroup <- as.factor(dataDGE$samples$subGroup)


###### Filtering lowly expressed reads

designDE <- model.matrix(~0 + subGroup, data = dataDGE$samples)
colnames(designDE) <- gsub(colnames(designDE), pattern = "subGroup", replacement = "")
keep <- filterByExpr.DGEList(dataDGE, designDE, min.count=10)
table(keep)
dataDGE <- dataDGE[keep,,keep.lib.sizes = F]
dataDGE <- calcNormFactors(dataDGE)
dataDGE <- estimateDisp(dataDGE, design = designDE)


######### Visualizing separation by MDS plot of counts #######

mds_BL6 <- plotMDS(dataDGE$counts, plot = FALSE,top = 500)

mds_df <- data.frame(X = mds_BL6$x, Y = mds_BL6$y,
                     diet = dataDGE$samples$Diet,
                     time = dataDGE$samples$Hours.from.lights.off,
                     photoperiod = dataDGE$samples$Photoperiod)

photoperiod_mds <- ggplot(mds_df, aes(X, Y)) + geom_text(aes(label=time, color=photoperiod), size=3) +
  theme_bw(base_size=10) +   scale_color_manual(values=Photocolour) + 
  theme(aspect.ratio=1) 

time_mds <- ggplot(mds_df, aes(X, Y)) + geom_point(aes(color=time), size=3) +
  theme_bw(base_size=10) + 
  theme(aspect.ratio=1) 

diet_mds <- ggplot(mds_df, aes(X, Y)) + geom_point(aes(color=diet), size=3) +
  theme_bw(base_size=10) + 
  theme(aspect.ratio=1) 

# separation by diet is most clear, no very obvious outliers

BL6_MDS_plots <- photoperiod_mds|time_mds|diet_mds
BL6_MDS_plots

# Design contrasts
contrastsDE <- makeContrasts(Diet = (SL_HFD + EL_HFD + LL_HFD)/3 - 
                               (SL_LFD + EL_LFD + LL_LFD)/3,
                             SLvsEL = (SL_HFD + SL_LFD)/2 - (EL_HFD + EL_LFD)/2,
                             LLvsEL = (LL_HFD + LL_LFD)/2 - (EL_HFD + EL_LFD)/2,
                             SLvsLL = (SL_HFD + SL_LFD)/2 - (LL_HFD + LL_LFD)/2,
                             SLvsELinHFD = SL_HFD - EL_HFD,
                             SLvsELinLFD = SL_LFD - EL_LFD,
                             LLvsELinHFD = LL_HFD - EL_HFD,
                             LLvsELinLFD = LL_LFD - EL_LFD,
                             SLvsLLinHFD = SL_HFD - LL_HFD,
                             SLvsLLinLFD = SL_LFD - LL_LFD,
                             levels = designDE)

# Voom and limma
voomDE <- voomWithQualityWeights(dataDGE, designDE, plot=TRUE)

fit <- lmFit(voomDE, designDE)

dataDGE[["DEfit"]] <- fit

resDE <- lapply(colnames(contrastsDE), function(cont){
  
  fitCont <- contrasts.fit(fit, contrasts = contrastsDE[,cont])
  fitCont <- eBayes(fitCont)
  topTable(fitCont, number = Inf, sort.by = "none")
  
})
names(resDE) <- colnames(contrastsDE)

######### 

# Add in ST time to metadata
dataDGE$samples <- dplyr::select(dataDGE$samples,!group)
dataDGE$samples <- dataDGE$samples %>% mutate(dataDGE$samples, STime = case_when(Photoperiod %in% "SL" ~ dataDGE$samples$Hours.from.lights.off-4.5,Photoperiod %in% "EL" ~ dataDGE$samples$Hours.from.lights.off-3,Photoperiod %in% "LL" ~ dataDGE$samples$Hours.from.lights.off-1.5)) %>% mutate(ST=sapply(STime,timer)) %>% dplyr::select(-STime)
dataDGE$samples <- merge(dataDGE$samples,Pooled.Metadata[c("Sample.ID","ExT")],by = "Sample.ID")

# Make log cpm table and tidy cpm table #
cpm_table <- cpm(dataDGE, log=TRUE)
cpm_tibble <- as_tibble(cpm_table, rownames = "genes")
tidy_cpm_table <- pivot_longer(cpm_tibble,"0101_1":"0101_123",names_to = "Sample.ID",values_to = "LogCPM")
tidy_cpm_table <- inner_join(dataDGE$samples,tidy_cpm_table,by="Sample.ID")
tidy_cpm_table <- inner_join(dataDGE$genes,tidy_cpm_table,by="genes")
tidy_cpm_table <- dplyr::select(tidy_cpm_table, !c("lib.size":"SCOP..","Project":"Cell.type","Sex":"Sequencing.Lane"))
tidy_cpm_table$Diet <- factor(tidy_cpm_table$Diet, levels = c("LFD", "HFD"))
tidy_cpm_table$Photoperiod <- factor(tidy_cpm_table$Photoperiod, levels = c("SL", "EL", "LL"))

### add average expression to dataDGE genes
identical(dataDGE$genes$genes,rownames(cpm_table))
dataDGE$genes <- mutate(dataDGE$genes,"average_expression"=rowMeans(cpm_table)) 

#Make summary tables for plotting 
summary_diet <- summarySE(tidy_cpm_table, measurevar="LogCPM", groupvars=c("genes","Diet","Photoperiod","zTime","Hours.from.lights.off","ExT","Symbol","ST"))
summary <- summarySE(tidy_cpm_table, measurevar="LogCPM", groupvars=c("genes","Photoperiod","zTime","Hours.from.lights.off","ExT","Symbol","ST"))

###########################################
# Rhythmicity Analysis with compareRhythms
##########################################

# Comparing between LFD and HFD for each photoperiod on Liver RNAseq
############################################################

# Compare Rhythms function for each diet individually 
DiffRhythm_diet <- function(pp){
  CR_meta_diet=subset(dataDGE$samples,Photoperiod==pp)
  CR_meta_diet=dplyr::rename(CR_meta_diet,group="Diet",time="zTime")
  CR_meta_diet=dplyr::select(CR_meta_diet, c("Sample.ID","group","time"))
  CR_meta_diet$group <- factor(CR_meta_diet$group, levels=c("LFD","HFD"))
  CRcounts <- subset(dataDGE$counts, select = CR_meta_diet$Sample.ID)
  CR_diet=compareRhythms::compareRhythms(CRcounts,CR_meta_diet,just_classify=FALSE,method =  "voom",amp_cutoff = 0.5,rhythm_fdr = 0.01,compare_fdr = 0.05, outliers=TRUE)
  CR_diet=CR_diet %>% add_column(LFD_Phase_h = (CR_diet$LFD_phase*24/(2*pi)))  
  CR_diet=CR_diet %>% add_column(HFD_Phase_h = (CR_diet$HFD_phase*24/(2*pi))) 
  CR_diet=CR_diet %>% add_column(amp_difference = (CR_diet$HFD_amp-CR_diet$LFD_amp))
  CR_diet=dplyr::select(CR_diet,!c("LFD_phase","HFD_phase"))
  CR_diet=mutate(CR_diet, phase_difference=NA)
  for (i in 1:length(CR_diet$LFD_Phase_h)){
    CR_diet$phase_difference[i] = hourCalculator(CR_diet$LFD_Phase_h[i],CR_diet$HFD_Phase_h[i])}
  return(CR_diet)}

# Using compareRhythms
SL_diet <- DiffRhythm_diet("SL")
EL_diet <- DiffRhythm_diet("EL")
LL_diet <- DiffRhythm_diet("LL")

table(SL_diet$category)
table(EL_diet$category)
table(LL_diet$category)

#Finding genes that are different in all photoperiods between LFD and HFD 
SL_diffm <- dplyr::select(SL_diet,c("id","category")) %>% dplyr::rename(category_SL=category)
EL_diffm <- dplyr::select(EL_diet, c("id","category")) %>% dplyr::rename(category_EL=category)
LL_diffm <- dplyr::select(LL_diet, c("id","category")) %>% dplyr::rename(category_LL=category)
diet_diffm <- SL_diffm %>% merge(EL_diffm,by="id",all = T) %>% merge(LL_diffm,by="id",all = T)

#Adding stats to table#
tidy_cpm_table_stats <- tidy_cpm_table %>% full_join(SL_diet[c("id","category")],by=c("genes"="id")) %>% dplyr::rename("category_SL"="category") %>% full_join(EL_diet[c("id","category")],by=c("genes"="id")) %>% dplyr::rename("category_EL"="category") %>% full_join(LL_diet[c("id","category")],by=c("genes"="id")) %>% dplyr::rename("category_LL"="category")

summary_diet_stats <-  summary_diet %>% full_join(SL_diet[c("id","category")],by=c("genes"="id")) %>% dplyr::rename("category_SL"="category") %>% full_join(EL_diet[c("id","category")],by=c("genes"="id")) %>% dplyr::rename("category_EL"="category") %>% full_join(LL_diet[c("id","category")],by=c("genes"="id")) %>% dplyr::rename("category_LL"="category")

# bar plots of diff rhythmic expression between LFD and HFD in liver in the different photoperiods
#############################################################################

#### bar plot of DODR significance 
tables_diet <- list(table(SL_diet$category),table(EL_diet$category),table(LL_diet$category))
c_names_diet <-  unique(unlist(SL_diet$category))
categories_diet <- do.call(rbind, lapply(tables_diet, `[`, c_names_diet)) %>% replace(is.na(.), 0)
rownames(categories_diet) <- c("SL","EL","LL") 
categories_diet <- as.data.frame(categories_diet) %>% rownames_to_column(var = "photoperiod")
categories_diet <- pivot_longer(as.data.frame(categories_diet),cols = 2:5,names_to = "category")
categories_diet$category <- factor(categories_diet$category, levels = c("change","gain","loss","same"))
categories_diet$photoperiod <- factor(categories_diet$photoperiod, levels = c("SL","EL","LL"))

#### Figure S4A ######

category_diet_bar_plot <- ggplot(categories_diet,aes(y=category,x=value,fill=photoperiod))+geom_col(position = "dodge",color="black")+scale_fill_manual(values=Photocolour)+theme_classic()+labs(y="",x="Number of Transcripts",fill="")+ theme(legend.position="top")+ geom_text(aes(label=value), position=position_dodge(width=0.9), hjust=-0.5)+ scale_x_continuous(trans='log10',limits = c(-1, 10000))

######## Phase and amplitude density of liver transcripts between LFD and HFD at different photo periods ###########
##############################################################################################

density_plots_amp_function <- function(pp_sheet,pp_color,pp){
  graph <- pp_sheet %>% filter(LFD_amp > 0 & HFD_amp > 0) %>% ggplot()+geom_density(aes(x=LFD_amp),size=1,color=pp_color)+geom_density(aes(x=HFD_amp),size=1,color=pp_color,linetype="dashed")+theme_minimal()+labs(x="Log2 Amplitude")+ggtitle(pp)+scale_x_continuous(limits = c(0,5))+scale_y_continuous(limits = c(0,1.5))+ theme(plot.title = element_text(hjust = 0.5))
  return(graph)}

density_plots_amp <- density_plots_amp_function(SL_diet,"#56B4E9","SL")/density_plots_amp_function(EL_diet,"black","EL")/density_plots_amp_function(LL_diet,"#FFCC33","LL")

density_plots_phase_function <- function(pp_sheet,pp_color,pp){
  graph <- pp_sheet %>% filter(LFD_amp > 1 & HFD_amp > 1) %>% ggplot(aes(x=-phase_difference))+geom_histogram(aes(y=..density..),alpha=0.5,fill=pp_color)+geom_density(size=1,color=pp_color)+geom_vline(xintercept = 0,color="black",linetype="dotted")+theme_minimal()+labs(x="phase diff LFD vs HFD")+scale_x_continuous(limits = c(-5,5))+scale_y_continuous(limits = c(0,0.5))+ggtitle(pp) + theme(plot.title = element_text(hjust = 0.5))
  return(graph)}

density_plots_phase <- density_plots_phase_function(SL_diet,"#56B4E9","SL")/density_plots_phase_function(EL_diet,"black","EL")/density_plots_phase_function(LL_diet,"#FFCC33","LL")

#### Figure S4B ######

density_plots_amp|density_plots_phase

# Comparing clock gene expression between LFD and HFD in liver (clock panel)
######################################################################
clock_panel_ensembl <- dplyr::filter(dataDGE$genes,dataDGE$genes$Symbol %in% c("Per2","Per1","Nr1d1","Npas2","Dbp","Arntl"))

filter(diet_diffm,diet_diffm$id %in% clock_panel_ensembl$genes) %>% inner_join(clock_panel_ensembl,by = c("id"="genes"))

summary_diet_stats$category_SL <- factor(summary_diet_stats$category_SL,levels = c("same","change","gain","loss"))
diet_clock_panel_liver <- function(pp,colour,title,category){
  graph <- summary_diet_stats %>% filter(Symbol %in% clock_panel) %>% filter(Photoperiod==pp) %>% 
    ggplot(aes_string(x="zTime",y="LogCPM",color=category))+geom_point(aes(shape=Diet),alpha=0.5)+geom_linerange(aes(ymin=LogCPM, ymax=LogCPM+se),alpha=0.5)+scale_color_manual(values=c(colour,"red"))+geom_smooth(aes(linetype=Diet),se=FALSE,span=0.5, method="lm", formula=y~cos(pi*x/12) + sin(pi*x/12))+theme_classic()+facet_wrap(~Symbol, scale="free_y",ncol = 3)+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,face = "bold"),strip.text.x = element_text(face = "italic")) + guides(shape="none",linetype="none",color="none")+scale_x_continuous(breaks=c(0,4,8,12,16,20))+labs(x=bquote(ZT[From_Lights_On]))
  return(graph)}

SL_diet_clock_plot <-  diet_clock_panel_liver("SL","#56B4E9","Short Light","category_SL")
EL_diet_clock_plot <-  diet_clock_panel_liver("EL","#000000","Equal Light","category_EL")
LL_diet_clock_plot <-  diet_clock_panel_liver("LL","#FFCC33","Long Light","category_LL")

######## Figure S4C ##########

lfd_vs_hfd_liver_clocks <- SL_diet_clock_plot/EL_diet_clock_plot/LL_diet_clock_plot

# Comparing clock gene expression between LFD and HFD in other tissues
#############################################################################

# load in Matrix of scaled qPCR data and metadata
qPCR_matrix <- read_csv("/Users/cnh943/Documents/Work/CBMR/Projects/Seasonal Study/RNA/R analysis/Seasonal qPCR/qPCR_matrix.csv")
qPCR_metadata <- read_delim("/Users/cnh943/Documents/Work/CBMR/Projects/Seasonal Study/RNA/R analysis/Seasonal qPCR/qPCR_metadata.csv",";", escape_double = FALSE, trim_ws = TRUE)

# Add in ST time to metadata
qPCR_metadata <-  qPCR_metadata %>% mutate(qPCR_metadata, STime = case_when(photoperiod %in% "SL" ~ qPCR_metadata$zt_off-4.5,photoperiod %in% "EL" ~ qPCR_metadata$zt_off-3,photoperiod %in% "LL" ~ qPCR_metadata$zt_off-1.5)) %>% mutate(ST=sapply(STime,timer)) %>% dplyr::select(-STime) %>% merge(y=Pooled.Metadata[c("Subject.Name","ExT")],by.x = "sample_id",by.y = "Subject.Name")

# Making qPCR matrix tidy
qPCR_tidy <- pivot_longer(qPCR_matrix,"1_1":"24_6",names_to = "sample_id",values_to = "LogFC")
qPCR_tidy <- inner_join(qPCR_metadata,qPCR_tidy,by="sample_id")
qPCR_tidy$photoperiod <- factor(qPCR_tidy$photoperiod, levels = c("SL", "EL", "LL"))
qPCR_tidy$diet <- factor(qPCR_tidy$diet, levels = c("LFD","HFD"))
qPCR_tidy$Tissue <- factor(qPCR_tidy$Tissue, levels = c("Hypo","Quad","WAT","BAT"))
qPCR_tidy$Gene <- factor(qPCR_tidy$Gene)
qPCR_tidy <- mutate(qPCR_tidy,tissue_gene = str_c(Tissue,"_",Gene))

qPCR_summary_diet <- summarySE(qPCR_tidy, measurevar="LogFC", groupvars=c("Gene","diet","photoperiod","zt_on","zt_off","ExT","Tissue","ST"))

#comparing between diets for each photoperiod with Compare Rhythms
DiffRhythm_qpcr=function(pp,meth){
  CRmatqpcr=subset(qPCR_metadata,photoperiod==pp)
  CRmatqpcr=dplyr::rename(CRmatqpcr,group=diet,time=zt_on)
  CRmatqpcr$group <- factor(CRmatqpcr$group,levels=c("LFD","HFD"))
  clock_genes=filter(qPCR_matrix, Gene %in% clock_panel )
  CRpcr=unite(clock_genes,col = Tissue_Gene, c("Tissue","Gene"),sep = "_")
  CRpcr= column_to_rownames(CRpcr, var="Tissue_Gene")
  CRpcr=subset(CRpcr, select = CRmatqpcr$sample_id)
  CRpcr=as.matrix(CRpcr)
  Rhythmpcr=compareRhythms::compareRhythms(CRpcr,CRmatqpcr,just_classify=FALSE,method = meth,amp_cutoff = 0.5,compare_fdr = 0.05,schwarz_wt_cutoff = 0.3,rhythm_fdr = 0.01)
  Rhythmpcr=Rhythmpcr %>% add_column(LFD_Phase_h = (Rhythmpcr$LFD_phase*24/(2*pi)))  
  Rhythmpcr=Rhythmpcr %>% add_column(HFD_Phase_h = (Rhythmpcr$HFD_phase*24/(2*pi))) 
  Rhythmpcr=Rhythmpcr %>% add_column(amp_difference = (Rhythmpcr$HFD_amp-Rhythmpcr$LFD_amp))
  Rhythmpcr=dplyr::select(Rhythmpcr,!c("LFD_phase","HFD_phase"))
  Rhythmpcr=mutate(Rhythmpcr, phase_difference=NA)
  for (i in 1:length(Rhythmpcr$LFD_Phase_h)){
    Rhythmpcr$phase_difference[i] = hourCalculator(Rhythmpcr$LFD_Phase_h[i],Rhythmpcr$HFD_Phase_h[i])}
  return(Rhythmpcr)}

SL_model_pcr <- DiffRhythm_qpcr("SL","mod_sel")
EL_model_pcr <- DiffRhythm_qpcr("EL","mod_sel")
LL_model_pcr <- DiffRhythm_qpcr("LL","mod_sel") 

SL_dodr_pcr <- DiffRhythm_qpcr("SL","dodr")

#Adding stats to table#

qPCR_tidy_stats <-  qPCR_tidy %>% unite("id", Tissue:Gene, sep = "_",remove = F) %>% full_join(SL_model_pcr[c("id","category")],by="id") %>% dplyr::rename("category_SL"="category") %>% full_join(EL_model_pcr[c("id","category")],by="id") %>% dplyr::rename("category_EL"="category") %>% full_join(LL_model_pcr[c("id","category")],by="id") %>% dplyr::rename("category_LL"="category")

qPCR_summary_diet <-  qPCR_summary_diet %>% unite("id", c("Tissue","Gene"), sep = "_",remove = F) %>% full_join(SL_model_pcr[c("id","category")],by="id") %>% dplyr::rename("category_SL"="category") %>% full_join(EL_model_pcr[c("id","category")],by="id") %>% dplyr::rename("category_EL"="category") %>% full_join(LL_model_pcr[c("id","category")],by="id") %>% dplyr::rename("category_LL"="category")

# Clock panel

diet_clock_panel_tissues <- function(pp,colour,title,category){
  graph <- qPCR_summary_diet %>% filter(Gene %in% clock_panel) %>% filter(photoperiod==pp) %>% 
    ggplot(aes_string(x="zt_on",y="LogFC",color=category))+geom_point(aes(shape=diet),alpha=0.5)+geom_linerange(aes(ymin=LogFC, ymax=LogFC+se),alpha=0.5)+scale_color_manual(values=c("gray",colour,"red"))+geom_smooth(aes(linetype=diet),se=FALSE,span=0.5, method="lm", formula=y~cos(pi*x/12) + sin(pi*x/12))+theme_classic()+facet_grid(Tissue~Gene, scale="free_y")+ggtitle(title)+theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"),strip.text.x = element_text(face = "italic")) + guides(shape="none",linetype="none")+scale_x_continuous(breaks=c(0,4,8,12,16,20))+labs(x=bquote(ZT[From_Lights_On]))
  return(graph)}

SL_diet_clock_plot_tissues <-  diet_clock_panel_tissues("SL","#56B4E9","Short Light","category_SL")
EL_diet_clock_plot_tissues <-  diet_clock_panel_tissues("EL","#000000","Equal Light","category_EL")
LL_diet_clock_plot_tissues <-  diet_clock_panel_tissues("LL","#FFCC33","Long Light","category_LL")

######## Figure S4D ##########

lfd_vs_hfd_tissue_clocks <- SL_diet_clock_plot_tissues/EL_diet_clock_plot_tissues/LL_diet_clock_plot_tissues

#### FIGURE S4 (patchwork) #######
######################

(category_diet_bar_plot/(density_plots_amp|density_plots_phase))|lfd_vs_hfd_liver_clocks|lfd_vs_hfd_tissue_clocks+plot_layout(widths=c(1,1,2))

##########################
####### FIGURE 3 #########
##########################

########### Clock gene panel liver
#############################################################################

clock_liver_panel <- function(t,xaxis){
  g <- summary %>% filter(Symbol %in% clock_panel) %>%
    ggplot(aes_string(x=t,y="LogCPM",color="Photoperiod"))+geom_linerange(aes(ymin=LogCPM, ymax=LogCPM+se))+geom_point()+geom_smooth(span=0.5,se=FALSE,method="lm", formula=y~cos(pi*x/12) + sin(pi*x/12))+scale_color_manual(values=Photocolour)+theme_classic()+facet_wrap(~Symbol, scale="free_y")+theme(plot.title = element_text(hjust = 0.5,face = "bold"),strip.text.x = element_text(face = "italic"),panel.border = element_rect(colour = "black", fill=NA, size=1))+scale_x_continuous(breaks=c(0,4,8,12,16,20))+ labs(x=xaxis)
  return(g)}

liver_on <- clock_liver_panel("zTime",bquote(ZT[From_Lights_On]))
liver_off <- clock_liver_panel("Hours.from.lights.off",bquote(ZT[From_Lights_Off]))
liver_ExT <- clock_liver_panel("ExT","ExT")
liver_ST <- clock_liver_panel("ST","ST")

######## Figure 3C ##########

liver_clock <- (liver_on/liver_off/liver_ExT/liver_ST)+plot_layout(guides="collect")

#### Clock gene panel all tissues
#############################################################################

#Finding phase and amplitude for other tissues 

DiffRhythm_qpcrSLvsEL <- function(ztime,meth){
  CRmatqpcr=subset(qPCR_metadata,photoperiod!="LL")
  CRmatqpcr=dplyr::rename(CRmatqpcr,group=photoperiod,time=ztime)
  CRmatqpcr$group <- factor(CRmatqpcr$group,levels=c("SL","EL"))
  clock_genes=filter(qPCR_matrix, Gene %in% clock_panel )
  CRpcr=unite(clock_genes,col = Tissue_Gene, c("Tissue","Gene"),sep = "_")
  CRpcr= column_to_rownames(CRpcr, var="Tissue_Gene")
  CRpcr=subset(CRpcr, select = CRmatqpcr$sample_id)
  CRpcr=as.matrix(CRpcr)
  Rhythmpcr=compareRhythms::compareRhythms(CRpcr,CRmatqpcr,just_classify=FALSE,method = meth,amp_cutoff = 0.5,compare_fdr = 0.05,schwarz_wt_cutoff = 0.3)
  Rhythmpcr=Rhythmpcr %>% add_column(SL_Phase_h = (Rhythmpcr$SL_phase*24/(2*pi)))  
  Rhythmpcr=Rhythmpcr %>% add_column(EL_Phase_h = (Rhythmpcr$EL_phase*24/(2*pi))) 
  Rhythmpcr=Rhythmpcr %>% add_column(amp_difference = (Rhythmpcr$SL_amp-Rhythmpcr$EL_amp))
  Rhythmpcr=dplyr::select(Rhythmpcr,!c("SL_phase","EL_phase"))
  Rhythmpcr=mutate(Rhythmpcr, phase_difference=NA)
  for (i in 1:length(Rhythmpcr$SL_Phase_h)){
    Rhythmpcr$phase_difference[i] = hourCalculator(Rhythmpcr$SL_Phase_h[i],Rhythmpcr$EL_Phase_h[i])}
  return(Rhythmpcr)}

qpcr_dodron_SLvsEL <- DiffRhythm_qpcrSLvsEL("zt_on","dodr") %>%  separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodroff_SLvsEL <- DiffRhythm_qpcrSLvsEL("zt_off","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodrExT_SLvsEL <- DiffRhythm_qpcrSLvsEL("ExT","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodrST_SLvsEL <- DiffRhythm_qpcrSLvsEL("ST","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")

DiffRhythm_qpcrLLvsEL <- function(ztime,meth){
  CRmatqpcr=subset(qPCR_metadata,photoperiod!="SL")
  CRmatqpcr=dplyr::rename(CRmatqpcr,group=photoperiod,time=ztime)
  CRmatqpcr$group <- factor(CRmatqpcr$group,levels=c("LL","EL"))
  clock_genes=filter(qPCR_matrix, Gene %in% clock_panel )
  CRpcr=unite(clock_genes,col = Tissue_Gene, c("Tissue","Gene"),sep = "_")
  CRpcr= column_to_rownames(CRpcr, var="Tissue_Gene")
  CRpcr=subset(CRpcr, select = CRmatqpcr$sample_id)
  CRpcr=as.matrix(CRpcr)
  Rhythmpcr=compareRhythms::compareRhythms(CRpcr,CRmatqpcr,just_classify=FALSE,method = meth,amp_cutoff = 0.5,compare_fdr = 0.05,schwarz_wt_cutoff = 0.3)
  Rhythmpcr=Rhythmpcr %>% add_column(LL_Phase_h = (Rhythmpcr$LL_phase*24/(2*pi)))  
  Rhythmpcr=Rhythmpcr %>% add_column(EL_Phase_h = (Rhythmpcr$EL_phase*24/(2*pi))) 
  Rhythmpcr=Rhythmpcr %>% add_column(amp_difference = (Rhythmpcr$LL_amp-Rhythmpcr$EL_amp))
  Rhythmpcr=dplyr::select(Rhythmpcr,!c("LL_phase","EL_phase"))
  Rhythmpcr=mutate(Rhythmpcr, phase_difference=NA)
  for (i in 1:length(Rhythmpcr$LL_Phase_h)){
    Rhythmpcr$phase_difference[i] = hourCalculator(Rhythmpcr$LL_Phase_h[i],Rhythmpcr$EL_Phase_h[i])}
  return(Rhythmpcr)}

qpcr_dodron_LLvsEL <- DiffRhythm_qpcrLLvsEL("zt_on","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodroff_LLvsEL <- DiffRhythm_qpcrLLvsEL("zt_off","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodrExT_LLvsEL <- DiffRhythm_qpcrLLvsEL("ExT","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodrST_LLvsEL <- DiffRhythm_qpcrLLvsEL("ST","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")

DiffRhythm_qpcrSLvsLL <- function(ztime,meth){
  CRmatqpcr=subset(qPCR_metadata,photoperiod!="EL")
  CRmatqpcr=dplyr::rename(CRmatqpcr,group=photoperiod,time=ztime)
  CRmatqpcr$group <- factor(CRmatqpcr$group,levels=c("SL","LL"))
  clock_genes=filter(qPCR_matrix, Gene %in% clock_panel )
  CRpcr=unite(clock_genes,col = Tissue_Gene, c("Tissue","Gene"),sep = "_")
  CRpcr= column_to_rownames(CRpcr, var="Tissue_Gene")
  CRpcr=subset(CRpcr, select = CRmatqpcr$sample_id)
  CRpcr=as.matrix(CRpcr)
  Rhythmpcr=compareRhythms::compareRhythms(CRpcr,CRmatqpcr,just_classify=FALSE,method = meth,amp_cutoff = 0.5,compare_fdr = 0.05,schwarz_wt_cutoff = 0.3)
  Rhythmpcr=Rhythmpcr %>% add_column(SL_Phase_h = (Rhythmpcr$SL_phase*24/(2*pi)))  
  Rhythmpcr=Rhythmpcr %>% add_column(LL_Phase_h = (Rhythmpcr$LL_phase*24/(2*pi))) 
  Rhythmpcr=Rhythmpcr %>% add_column(amp_difference = (Rhythmpcr$SL_amp-Rhythmpcr$LL_amp))
  Rhythmpcr=dplyr::select(Rhythmpcr,!c("SL_phase","LL_phase"))
  Rhythmpcr=mutate(Rhythmpcr, phase_difference=NA)
  for (i in 1:length(Rhythmpcr$SL_Phase_h)){
    Rhythmpcr$phase_difference[i] = hourCalculator(Rhythmpcr$SL_Phase_h[i],Rhythmpcr$LL_Phase_h[i])}
  return(Rhythmpcr)}

qpcr_dodron_SLvsLL <- DiffRhythm_qpcrSLvsLL("zt_on","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodroff_SLvsLL <- DiffRhythm_qpcrSLvsLL("zt_off","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodrExT_SLvsLL <- DiffRhythm_qpcrSLvsLL("ExT","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")
qpcr_dodrST_SLvsLL <- DiffRhythm_qpcrSLvsLL("ST","dodr") %>% separate(col = "id",into = c("tissue","gene"),sep = "_")

##### Testing differential rhythmicity in liver RNAseq using compare rhythms between photoperiods combining diets ##########
############################################################################################################

#voom-dodr of liver data

DiffRhythm <- function(ZT,pp1,pp2){
  CRmatrix=subset(dataDGE$samples,Photoperiod==pp1|Photoperiod==pp2)
  CRmatrix=dplyr::rename(CRmatrix,group="Photoperiod",time=ZT,batch="Diet")
  CRmatrix=dplyr::select(CRmatrix, c("Sample.ID","group","time","batch"))
  CRmatrix$group <- as.factor(CRmatrix$group)
  CRmatrix$group <- droplevels(CRmatrix$group)
  CRmatrix$group <-relevel(CRmatrix$group, pp2)
  CRmatrix$batch <- as.factor(CRmatrix$batch)
  CRcounts <- subset(dataDGE$counts, select = CRmatrix$Sample.ID)
  Rhythm=compareRhythms::compareRhythms(CRcounts,CRmatrix,just_classify=FALSE,method =  "voom",amp_cutoff = 0.5,rhythm_fdr = 0.01,compare_fdr = 0.05, outliers=TRUE)
  return(Rhythm)}
addamp <- function(sheet,pp1,pp2){
  names(sheet) <- gsub(pp1,"p1",names(sheet))
  names(sheet) <- gsub(pp2,"p2",names(sheet))
  sheet=sheet %>% add_column(p1_Phase_h = (sheet$p1_phase*24/(2*pi)))
  sheet=sheet %>% add_column(p2_Phase_h = (sheet$p2_phase*24/(2*pi)))
  sheet=sheet %>% add_column(amp_difference = (sheet$p1_amp-sheet$p2_amp))
  sheet=dplyr::select(sheet,!c("p1_phase","p2_phase"))
  sheet=mutate(sheet, phase_difference=NA)
  for (i in 1:length(sheet$p1_Phase_h)){
    sheet$phase_difference[i] = hourCalculator(sheet$p1_Phase_h[i],sheet$p2_Phase_h[i])}
  names(sheet) <- gsub("p1",pp1,names(sheet))
  names(sheet) <- gsub("p2",pp2,names(sheet))
  return(sheet)}

######ZT Lights On#######

dodr_on_SLvsLL <- DiffRhythm('zTime',"SL","LL") %>% addamp("SL","LL")
dodr_on_SLvsEL <- DiffRhythm('zTime',"SL","EL") %>% addamp("SL","EL")
dodr_on_LLvsEL <- DiffRhythm('zTime',"LL","EL") %>% addamp("LL","EL")

######ZT Lights Off#######

dodr_off_SLvsLL <- DiffRhythm("Hours.from.lights.off","SL","LL") %>% addamp("SL","LL")
dodr_off_SLvsEL <- DiffRhythm("Hours.from.lights.off","SL","EL") %>% addamp("SL","EL")
dodr_off_LLvsEL <- DiffRhythm("Hours.from.lights.off","LL","EL") %>% addamp("LL","EL")

##########ExT#########

dodr_ExT_SLvsLL <- DiffRhythm('ExT',"SL","LL") %>% addamp("SL","LL")
dodr_ExT_SLvsEL <- DiffRhythm('ExT',"SL","EL") %>% addamp("SL","EL")
dodr_ExT_LLvsEL <- DiffRhythm('ExT',"LL","EL") %>% addamp("LL","EL")

###### ST #######

dodr_ST_SLvsLL <- DiffRhythm('ST',"SL","LL") %>% addamp("SL","LL")
dodr_ST_SLvsEL <- DiffRhythm('ST',"SL","EL") %>% addamp("SL","EL")
dodr_ST_LLvsEL <- DiffRhythm('ST',"LL","EL") %>% addamp("LL","EL")

table(dodr_ST_SLvsLL$category)
table(dodr_ST_SLvsEL$category)
table(dodr_ST_LLvsEL$category)

#Adding liver data to qpcr data

phase_converter <- function(liver_sheet,qpcr_sheet,pp){
  a1=filter(liver_sheet,liver_sheet$id %in% clock_panel_ensembl$genes)
  a1= left_join(a1,dataDGE$genes,by = c("id" = "genes")) %>% mutate(tissue="Liver") %>% dplyr::select(c("tissue","Symbol","phase_difference","amp_difference")) %>% dplyr::rename("gene"= Symbol) 
  a2=dplyr::select(qpcr_sheet, c("tissue","gene","phase_difference","amp_difference"))
  a2=rbind(a2,a1)
  a2=mutate(a2,photoperiod=pp)
  return(a2)}

phasediff_SL_on <- phase_converter(dodr_on_SLvsEL,qpcr_dodron_SLvsEL,"SL")
phasediff_LL_on <- phase_converter(dodr_on_LLvsEL,qpcr_dodron_LLvsEL,"LL")
phasediff_SL_off <- phase_converter(dodr_off_SLvsEL,qpcr_dodroff_SLvsEL,"SL")
phasediff_LL_off <- phase_converter(dodr_off_LLvsEL,qpcr_dodroff_LLvsEL,"LL")
phasediff_SL_ExT <- phase_converter(dodr_ExT_SLvsEL,qpcr_dodrExT_SLvsEL,"SL")
phasediff_LL_ExT <- phase_converter(dodr_ExT_LLvsEL,qpcr_dodrExT_LLvsEL,"LL")
phasediff_SL_ST <- phase_converter(dodr_ST_SLvsEL,qpcr_dodrST_SLvsEL,"SL")
phasediff_LL_ST <- phase_converter(dodr_ST_LLvsEL,qpcr_dodrST_LLvsEL,"LL")

phasediff_on <- rbind(phasediff_SL_on,phasediff_LL_on)
phasediff_on$photoperiod <- factor(phasediff_on$photoperiod,levels = c("SL", "LL"))

phasediff_off <- rbind(phasediff_SL_off,phasediff_LL_off)
phasediff_off$photoperiod <- factor(phasediff_off$photoperiod,levels = c("SL", "LL"))

phasediff_ExT <- rbind(phasediff_SL_ExT,phasediff_LL_ExT)
phasediff_ExT$photoperiod <- factor(phasediff_ExT$photoperiod,levels = c("SL", "LL"))

phasediff_ST <- rbind(phasediff_SL_ST,phasediff_LL_ST)
phasediff_ST$photoperiod <- factor(phasediff_ST$photoperiod,levels = c("SL", "LL"))

tissue_phase_diff_plotter <- function(sheet,title){
  graph <- ggplot(sheet, aes(color=photoperiod,shape=tissue,x = phase_difference,y = gene))+geom_point(position = position_dodge(width = 0.30))+scale_color_manual(values = c("#56B4E9","#FFCC33"))+scale_x_continuous(limits = c(-6,6))+ labs(x = "Phase Diff from EL (h)", y = "") + ggtitle(title)+theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(face = "italic"))+geom_vline(xintercept = median(filter(sheet,photoperiod == "SL")$phase_difference),linetype="dashed",color="dark blue")+geom_vline(xintercept = median(filter(sheet,photoperiod == "LL")$phase_difference),linetype="dashed",color="#AA8E00")+guides(color="none")
return(graph)
  }

## Figure 3D ###

tissue_phase_difference <- tissue_phase_diff_plotter(phasediff_on,bquote(ZT[lights_on]))/tissue_phase_diff_plotter(phasediff_off,bquote(ZT[lights_off]))/tissue_phase_diff_plotter(phasediff_ExT,"ExT")/tissue_phase_diff_plotter(phasediff_ST,"ST")+ plot_layout(guides = "collect")

# Polar frequency of phase difference plot Liver RNA
#############################################################################

polar_histogram=function(df1,df2,title){
  a <- filter(df1,SL_amp>1&EL_amp>1)
  b <- filter(df2,LL_amp>1&EL_amp>1)
  graph <- ggplot() + geom_hline(yintercept = seq(0, 100, by = 25), colour = "grey90", size = 0.4
  ) + geom_vline(xintercept = seq(-9, 12, by = 3), colour = "grey90", size = 0.4
  ) + geom_histogram(
    bins = 80,
    fill = "#56B4E9",
    color = "black",
    data = a,
    aes(x = phase_difference)
  ) + geom_histogram(
    bins = 80,
    fill = "#FFCC33",
    color = "black",
    data = b,
    aes(x = phase_difference)
  ) + coord_polar(theta = "x", offset(pi)) + scale_x_continuous(limits = c(-12, 12),
                                                                breaks = c(0, 3, 6, 9, -3, -6, -9, 12)) + scale_y_continuous(limits = c(0, 100), breaks=c(25,50,75,100)
                                                                ) + theme(
                                                                  axis.text.y = element_blank(),
                                                                  axis.ticks.y = element_line(size = 0),
                                                                  axis.text.x = element_text(size = 15,face ="bold"),
                                                                  plot.title = element_text(hjust = 0.5),
                                                                  panel.grid.major = element_blank(),
                                                                  panel.background = element_blank()
                                                                ) + labs(x = "Phase Diff from EL (h)", y = "") + ggtitle(title)
  return(graph)}

polar_ZT <- polar_histogram(dodr_on_SLvsEL,dodr_on_LLvsEL,bquote(ZT[Lights_On]))
polar_ZToff <- polar_histogram(dodr_off_SLvsEL,dodr_off_LLvsEL,bquote(ZT[Lights_Off]))
polar_ExT <- polar_histogram(dodr_ExT_SLvsEL,dodr_ExT_LLvsEL,bquote("ExT"))
polar_ST <- polar_histogram(dodr_ST_SLvsEL,dodr_ST_LLvsEL,bquote("ST"))

## Figure 3B ###

polars <- polar_ZT/polar_ZToff/polar_ExT/polar_ST

######## FIGURE 3 ##########
########################################################

(((polar_ZT|liver_on)/(polar_ZToff|liver_off)/(polar_ExT|liver_ExT)/(polar_ST|liver_ST))|tissue_phase_difference) + plot_layout(widths = c(3,1))

(polars|liver_clock|tissue_phase_difference) + plot_layout(widths = c(1,2,1))

######## Finding genes with differential rhythmicity between SL and LL ########
############################################################################

SLvsLL_CR <- merge(dodr_ST_SLvsLL,dataDGE$genes,by.x="id",by.y="genes")

SLvsLL_same <- filter(SLvsLL_CR,category == "same")
SLvsLL_gain <- filter(SLvsLL_CR,category == "gain")
SLvsLL_loss <- filter(SLvsLL_CR,category == "loss")
SLvsLL_change <- filter(SLvsLL_CR,category == "change")

SLvsLL_DR <- filter(SLvsLL_CR,SLvsLL_CR$diff_rhythmic == T)
count(SLvsLL_DR$amp_difference > 0.5)
count(SLvsLL_DR$amp_difference < -0.5)

#### barplot of DODR significance 
tables=list(table(dodr_ST_SLvsLL$category),table(dodr_ST_SLvsEL$category),table(dodr_ST_LLvsEL$category))
c_names = unique(unlist(dodr_ST_SLvsLL$category))
categories <- do.call(rbind, lapply(tables, `[`, c_names)) %>% replace(is.na(.), 0)
rownames(categories) <- c("SLvsLL","SLvsEL","LLvsEL")
categories <- as.data.frame(categories) %>% rownames_to_column(var = "photoperiods")
categories <- pivot_longer(as.data.frame(categories),cols = 2:5,names_to = "category")
categories$category <- factor(categories$category, levels = c("change","gain","loss","same"))
categories$photoperiods <- factor(categories$photoperiods, levels = c("SLvsLL","SLvsEL","LLvsEL"))

##### Figure 4A ########
category_bar_plot <- ggplot(categories,aes(y=category,x=value,fill=photoperiods))+geom_col(position = "dodge",color="black")+scale_fill_manual(values=c("#9DC1A3","#56B4E9","#FFCC33"))+theme_classic()+labs(y="",x="Number of Transcripts",fill="")+ theme(legend.position="top")+ geom_text(aes(label=value), position=position_dodge(width=0.9), hjust=-0.5) + scale_x_continuous(trans='log10',limits = c(-1, 10000))

##### Density plot combined diets ########

density_amp_diets <- merge(dodr_ST_SLvsEL,dodr_ST_LLvsEL,by="id", all=T) %>% dplyr::select(c("id",SL_amp,EL_amp.x,LL_amp)) %>% dplyr::rename(EL_amp = EL_amp.x) %>% pivot_longer(cols = 2:4,values_to = "amp",names_to = "photoperiod")
density_amp_diets$photoperiod <- factor(density_amp_diets$photoperiod,levels = c("SL_amp","EL_amp","LL_amp"))

##### Figure 4B ########
density_plot_diet <- ggplot(density_amp_diets,aes(x=amp,color=photoperiod))+geom_density(size=1)+scale_color_manual(labels=c("SL","EL","LL"),values=Photocolour)+scale_x_continuous(limits = c(0, 5))+ labs(x="log2(amplitude)")+ theme(plot.title = element_text(hjust = 0.5))+ labs(fill = "Photoperiod",color="Photoperiod")+theme_classic()

#####Plotting amp vs expression

amp_expression <- SLvsLL_CR %>% filter(average_expression > 2) %>%  mutate(cutoff=ifelse(amp_difference < -0.5,"LL",(ifelse(amp_difference > 0.5,"SL","NO")))) %>%  mutate(label = ifelse(amp_difference > 1.3|amp_difference < -1.3,"YES","NO"))
amp_expression$cutoff <- factor(amp_expression$cutoff,levels = c("LL","NO","SL"))
amp_expression$category <- factor(amp_expression$category,levels = c("same","gain","loss","change"))
amp_expression$label <- factor(amp_expression$label,levels = c("YES","NO"))

##### Figure 4C ########
amp_vs_express <- ggplot(amp_expression, aes(x=amp_difference,y=average_expression,color=cutoff))+ geom_point(alpha=0.5)+ labs(x="log2 amplitude difference (SLvsLL) ",y="LogCPM",color="")+scale_color_manual(labels = c("Higher in LL", "Similar Amp", "Higher in SL"),values=c("#FFCC33","#000000", "#56B4E9"))+theme_classic()+theme(legend.position="top")+ geom_text( 
  data=amp_expression %>% filter(amp_expression$label=="YES"),
  aes(label=Symbol),nudge_x = 0.3,check_overlap = T,show.legend = F)

######## gene over representation analysis for genes different in amplitude

amp_loss <- filter(amp_expression,amp_difference < -0.5)
amp_loss_BP_GSE <- clusterProfiler::enrichGO(amp_loss$entrez,OrgDb = "org.Mm.eg.db",universe=amp_expression$entrez,ont="BP",pvalueCutoff = 0.1)
amp_gain <- filter(amp_expression,amp_difference > 0.5)
amp_gain_BP_GSE <- clusterProfiler::enrichGO(amp_gain$entrez,OrgDb = "org.Mm.eg.db",universe=amp_expression$entrez,ont="BP",pvalueCutoff = 0.1)

SLvsLL_change_BP <- enrichGO(gene = SLvsLL_change$entrez,OrgDb = "org.Mm.eg.db",universe=SLvsLL_CR$entrez,ont="BP",pvalueCutoff = 0.1)

##### Figure 4D ########
Merged_dotplot_BP <- merge_result(list(amp_loss = amp_loss_BP_GSE, rhythm_change = SLvsLL_change_BP, amp_gain = amp_gain_BP_GSE)) %>%
  dotplot(., showCategory=5)

#### FIGURE 4 (patchwork) #####

((category_bar_plot/density_plot_diet/amp_vs_express)|Merged_dotplot_BP)+plot_layout(widths = c(1,2))

#### Plotting genes with top amp difference

top_amp_in_SL <- SLvsLL_DR %>% filter(average_expression > 2) %>% slice_max(amp_difference,n=4)
top_amp_in_LL <- SLvsLL_DR %>% filter(average_expression > 2) %>% slice_min(amp_difference,n=4)
top_changed <- SLvsLL_DR %>% filter(average_expression > 2) %>% filter(category == "change") %>% slice_min(adj_p_val_DR,n=4)%>% slice_min(adj_p_val_LL_or_SL,n=4)

##### Figure 4E ########

(circadian_plotter_BL6(gene = genes,sheet = top_amp_in_LL$id,time = "ST",columns = 1,smoothing = cosinor_curve)|circadian_plotter_BL6(gene = genes,sheet = top_changed$id,time = "ST",columns = 1,smoothing = cosinor_curve)|circadian_plotter_BL6(gene = genes,sheet = top_amp_in_SL$id,time = "ST",columns = 1,smoothing = cosinor_curve))+plot_layout(guides = "collect")

#############################
###### C3H dataset ##########
#############################

#read in data
C3H_metadata <- read.csv("~/Documents/Work/CBMR/Projects/Seasonal Study/RNA/C3H liver RNA-seq/C3H Liver RNAseq/C3H_Metadata_Sheet.csv")
load("~/Documents/Work/CBMR/Projects/Seasonal Study/RNA/C3H liver RNA-seq/C3H Liver RNAseq/C3H.data.raw.Rdata")
C3H_counts <- C3H.data.raw
colnames(C3H_counts)=c("0143_10","0143_11","0143_12","0143_13","0143_14","0143_15","0143_16","0143_17","0143_18","0143_19","0143_01" ,"0143_20","0143_21","0143_22","0143_23","0143_24","0143_25","0143_26","0143_27","0143_28","0143_29","0143_02","0143_30","0143_31","0143_32","0143_03","0143_04","0143_05", "0143_06","0143_07","0143_08","0143_09")

####### Design #######
col_order <- C3H_metadata$sample_id
C3H_counts <- C3H_counts[, col_order]
C3H_DGE <- DGEList(counts = C3H_counts,samples =  C3H_metadata)
ztime <- factor(C3H_DGE$samples$ZToff)
photoperiod <- factor(C3H_DGE$samples$Photoperiod,levels = c("SL","LL"))
require(org.Mm.eg.db)
Symbol <- mapIds(org.Mm.eg.db, keys=rownames(C3H_DGE$counts), keytype="ENSEMBL",column = "SYMBOL")
ENTREZ <- mapIds(org.Mm.eg.db, keys=rownames(C3H_DGE$counts), keytype="ENSEMBL",column = "ENTREZID")
C3H_DGE$genes <- data.frame(ENSEMBL=rownames(C3H_DGE$counts),Symbol=Symbol,ENTREZ=ENTREZ)
C3H_design <- model.matrix(~0 + Photoperiod,data = C3H_DGE$samples)
colnames(C3H_design) <- gsub(colnames(C3H_design), pattern = "Photoperiod", replacement = "")
C3Hcontrasts <- makeContrasts(photo = SL - LL, levels = C3H_design)

##### Filter by expression ######

C3H_keep <- filterByExpr(C3H_DGE,design = C3H_design)
table(C3H_keep)
C3H_DGE <- C3H_DGE[C3H_keep,,keep.lib.sizes = F]
C3H_DGE <- calcNormFactors(C3H_DGE)
C3H_DGE <- estimateDisp(C3H_DGE, design = C3H_design)

# Voom and limma
voomDE_C3H <- voomWithQualityWeights(C3H_DGE, C3H_design, plot=TRUE)

fit_C3H <- lmFit(voomDE_C3H, C3H_design)

C3H_DGE[["DEfit"]] <- fit_C3H

C3HresDE <- lapply(colnames(C3Hcontrasts), function(cont){
  
  fitCont <- contrasts.fit(fit_C3H, contrasts = C3Hcontrasts[,cont])
  fitCont <- eBayes(fitCont)
  topTable(fitCont, number = Inf, sort.by = "none")
  
})
names(C3HresDE) <- colnames(C3Hcontrasts)

hist(C3HresDE$photo$P.Value)
C3H_cpm <- cpm(C3H_DGE, log=TRUE)

### MDS ###
C3H_mds <- plotMDS(C3H_DGE, top = 500, plot = F)

C3H_mds_df <- data.frame(X = C3H_mds$x, Y = C3H_mds$y,
                     time = C3H_DGE$samples$ZToff,
                     photoperiod = C3H_DGE$samples$Photoperiod)

ggplot(C3H_mds_df, aes(X, Y,color=time,shape=photoperiod)) + geom_point(size=3) +
  theme_bw(base_size=10) + 
  theme(aspect.ratio=1) 

# Add in ST and ExT time to metadata
C3H_DGE$samples <- C3H_DGE$samples %>% mutate(C3H_DGE$samples, STime = case_when(Photoperiod %in% "SL" ~ C3H_DGE$samples$ZToff-4.5,Photoperiod %in% "LL" ~ C3H_DGE$samples$ZToff-1.5)) %>% mutate(ST=sapply(STime,timer)) %>% dplyr::select(-STime)
C3H_DGE$samples <- C3H_DGE$samples %>% mutate(C3H_DGE$samples, ExTime = case_when(Photoperiod %in% "SL" ~ C3H_DGE$samples$ZToff-9,Photoperiod %in% "LL" ~ C3H_DGE$samples$ZToff-3)) %>% mutate(ExT=sapply(ExTime,timer)) %>% dplyr::select(-ExTime)

# Add in average expression to C3H_DGE$genes
identical(C3H_DGE$genes$ENSEMBL,rownames(C3H_cpm))
C3H_DGE$genes <- mutate(C3H_DGE$genes,"average_expression"=rowMeans(C3H_cpm)) 

#Making tidy dataset
#matrix to tibble 
library(tidybulk)
C3H_tidy <- as.data.frame(C3H_cpm)
C3H_tidy <- C3H_tidy %>% rownames_to_column() %>% pivot_longer(cols = c("0143_01":"0143_32"), names_to = "sample_id",values_to="LogCPM") %>% dplyr::rename("ENSEMBL"="rowname") %>% inner_join(C3H_DGE$samples,by = "sample_id") %>% inner_join(C3H_DGE$genes,by = "ENSEMBL")
C3H_tidy$Photoperiod <- factor(C3H_tidy$Photoperiod, levels = c("SL","LL"))
summary_C3H <- summarySE(C3H_tidy, measurevar="LogCPM", groupvars=c("ENSEMBL","Photoperiod","ZT","ZToff","Symbol","ST","ExT"))

#plotting
C3H_colors <- c("#915AF3","#F3A62B")

#compare Rhythms
C3H_DiffRhythm <- function(t){
  C3H_CRmeta <- dplyr::select(C3H_DGE$samples,-"group") %>% dplyr::rename(group="Photoperiod",time = all_of(t) ) %>% dplyr::select(c("sample_id","group","time"))
  C3H_CRmeta$group <- as.factor(C3H_CRmeta$group)
  C3H_CRmeta$group <- relevel(C3H_CRmeta$group,"SL","LL")
  Rhythm <- compareRhythms::compareRhythms(C3H_DGE$counts,C3H_CRmeta,just_classify=FALSE,method =  "voom",outliers = T,amp_cutoff = 0.5,rhythm_fdr = 0.01,compare_fdr = 0.05)
  Rhythm <- Rhythm %>% add_column(SL_Phase_h = (Rhythm$SL_phase*24/(2*pi)))  
  Rhythm <- Rhythm %>% add_column(LL_Phase_h = (Rhythm$LL_phase*24/(2*pi))) 
  Rhythm <- Rhythm %>% add_column(amp_difference = (Rhythm$SL_amp-Rhythm$LL_amp))
  Rhythm <- dplyr::select(Rhythm,-c("SL_phase","LL_phase"))
  Rhythm <- mutate(Rhythm, phase_difference=NA)
  for (i in 1:length(Rhythm$SL_Phase_h)){
    Rhythm$phase_difference[i] = hourCalculator(Rhythm$SL_Phase_h[i],Rhythm$LL_Phase_h[i])}
  return(Rhythm)}

C3H_on_dodr <- C3H_DiffRhythm("ZT")
C3H_off_dodr <- C3H_DiffRhythm("ZToff")
C3H_ExT_dodr <- C3H_DiffRhythm("ExT")
C3H_ST_dodr <- C3H_DiffRhythm("ST")

##### Figure S7B (Polar histograms) ##########
polar_histogram_C3H <- function(sheet){
  a <- filter(sheet,SL_amp>1&LL_amp>1)
  graph=ggplot() + geom_hline(yintercept = seq(0, 100, by = 25), colour = "grey90", size = 0.4
  ) + geom_vline(xintercept = seq(-9, 12, by = 3), colour = "grey90", size = 0.4
  ) + geom_histogram(
    bins = 60,
    fill = "grey",
    color = "black",
    data = a,
    aes(x = phase_difference)
  ) + coord_polar(theta = "x", offset(pi)) + scale_x_continuous(limits = c(-12, 12),
                                                                breaks = c(0, 3, 6, 9, -3, -6, -9, 12)) + scale_y_continuous(limits = c(0, 100), breaks=c(25,50,75,100)
                                                                ) + theme(
                                                                  axis.text.y = element_blank(),
                                                                  axis.ticks.y = element_line(size = 0),
                                                                  axis.text.x = element_text(size = 15,face ="bold"),
                                                                  plot.title = element_text(hjust = 0.5),
                                                                  panel.grid.major = element_blank(),
                                                                  panel.background = element_blank()
                                                                ) + labs(x = "Phase Diff SL vs LL", y = "")
  return(graph)}
C3H_polars <- polar_histogram_C3H(C3H_on_dodr)/polar_histogram_C3H(C3H_off_dodr)/polar_histogram_C3H(C3H_ExT_dodr)/polar_histogram_C3H(C3H_ST_dodr)

##### Figure S7C (circadian clock genes in liver) ##########
C3H_circplot_errorbars <- function(gene,sheet,time,columns){
  d <- filter(summary_C3H,{{gene}} %in% sheet)
  graph <- ggplot(d,aes_string(x = time, y = "LogCPM", color = "Photoperiod")) + geom_linerange(aes(ymin=LogCPM, ymax=LogCPM+se)) +geom_point()+ geom_smooth(
    se = FALSE,
    span = 0.5,
    method = "lm",
    formula = y ~ cos(pi * x / 12) + sin(pi * x / 12)
  ) + theme_classic() + facet_wrap( ~ Symbol, scale = "free_y",ncol = columns)+ scale_color_manual(values=C3H_colors)+scale_x_continuous(breaks=c(0,6,12,18,24),limits = c(0,24))+theme(strip.text = element_text(face = "italic"))+labs(y="LogCPM")
  return(graph)}
C3H_clock_sync <- C3H_circplot_errorbars(Symbol,clock_panel,"ZT",3)/C3H_circplot_errorbars(Symbol,clock_panel,"ZToff",3)/C3H_circplot_errorbars(Symbol,clock_panel,"ExT",3)/C3H_circplot_errorbars(Symbol,clock_panel,"ST",3)

##### Figure S7G (dodr barplot) ##########
table(C3H_ST_dodr$category)
C3H_dodr_table <- data.frame(category=c("change","gain","loss","same"),value=c(0,0,0,1885))
C3H_dodr_table$category <- factor(C3H_dodr_table$category, levels = c("change","gain","loss","same"))
C3H_bar_plot <- ggplot(C3H_dodr_table,aes(y=category,x=value))+geom_col(position = "dodge",color="black",fill="#9DC1A3")+theme_classic()+labs(y="",x="Number of Transcripts",fill="")+ theme(legend.position="top")+ geom_text(aes(label=value), position=position_dodge(width=0.9), hjust=-0.5) + scale_x_continuous(limits = c(0, 3000))

###### Figure S7H (Amp Density) ##########
density_ampSLvsLL_C3H <- dplyr::select(filter(C3H_ST_dodr,SL_amp !=0 & LL_amp !=0), c("SL_amp","LL_amp")) %>% reshape2::melt()
C3H_amp_density <- ggplot(density_ampSLvsLL_C3H,aes(x=value,color=variable))+geom_line(stat="density", alpha=0.7,size=1)+scale_color_manual(values=C3H_colors,labels=c("SL","LL"))+labs(x="log2(amplitude)")+theme_classic()+theme(legend.position = c(0.8, 0.8))+ labs(fill = "Photoperiod",color="Photoperiod")+scale_x_continuous(limits = c(0, 5))

##### Figure S7I (amp vs expression) #######

amp_expression_C3H <- inner_join(C3H_ST_dodr,C3H_DGE$genes[c("ENSEMBL","average_expression","Symbol","ENTREZ")],by=c("id"="ENSEMBL")) %>% filter(average_expression > 2) %>%  mutate(cutoff=ifelse(amp_difference < -0.5,"LL",(ifelse(amp_difference > 0.5,"SL","NO")))) %>%  mutate(label = ifelse(amp_difference > 1|amp_difference < -1,"YES","NO"))
amp_expression_C3H$cutoff <- factor(amp_expression_C3H$cutoff,levels = c("LL","NO","SL"))
C3H_amp_expression <- ggplot(amp_expression_C3H, aes(x=amp_difference,y=average_expression,color=cutoff))+ geom_point(alpha=0.5)+ labs(x="log2 amplitude difference (SLvsLL) ",y="LogCPM",color="")+scale_color_manual(labels = c("Higher in LL", "Similar Amp", "Higher in SL"),values=c("#F3A62B","#000000", "#915aF3"))+theme_classic()+theme(legend.position="top") + geom_text( 
  data=amp_expression_C3H %>% filter(label=="YES"),
  aes(label=Symbol),nudge_x = 0.3,check_overlap = T,show.legend = F)

#patchwork
C3H_polars|C3H_clock_sync

C3H_bar_plot|C3H_amp_density|C3H_amp_expression

###############################################
#### Figure S5 hypothalamus neuropeptides #####

neuro_genes <- c("Dio2","Npy","Agrp","Pomc","Pmch","Hcrt","Vip","Nms")

# determining outliers

(qPCR_summary_diet %>% filter(Gene %in% neuro_genes)%>% ggplot(aes(x=Gene,y=LogFC,color=photoperiod))+geom_boxplot()+scale_color_manual(values=Photocolour))/
(qPCR_summary_diet %>% filter(Gene %in% neuro_genes)%>% ggplot(aes(x=Gene,y=LogFC))+geom_boxplot())

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

qPCR_tidy <- qPCR_tidy %>% 
  group_by(Gene,Tissue) %>% 
  mutate("outlier" = is_outlier(LogFC))

pcr_outliers <- filter(qPCR_tidy, outlier == T)

qPCR_tidy_out <- qPCR_tidy %>% filter(outlier == F)

qPCR_summary_diet_out <- qPCR_tidy_out %>% summarySE(measurevar="LogFC", groupvars=c("Gene","diet","photoperiod","zt_on","zt_off","ExT","Tissue","ST"))

# testing rhythmicity in neuropeptides

DiffRhythm_neuroSLvsLL=function(ztime,meth){
  CRmatqpcr=subset(qPCR_metadata,photoperiod!="EL")
  CRmatqpcr=dplyr::rename(CRmatqpcr,group=photoperiod,time=ztime)
  CRmatqpcr$group <- factor(CRmatqpcr$group,levels=c("SL","LL"))
  n_genes=filter(qPCR_matrix, Gene %in% neuro_genes )
  CRpcr=unite(n_genes,col = Tissue_Gene, c("Tissue","Gene"),sep = "_")
  CRpcr= column_to_rownames(CRpcr, var="Tissue_Gene")
  CRpcr=subset(CRpcr, select = CRmatqpcr$sample_id)
  CRpcr=as.matrix(CRpcr)
  Rhythmpcr=compareRhythms::compareRhythms(CRpcr,CRmatqpcr,just_classify=FALSE,method = meth,amp_cutoff = 0.5,compare_fdr = 0.05,schwarz_wt_cutoff = 0.3)
  Rhythmpcr=Rhythmpcr %>% add_column(SL_Phase_h = (Rhythmpcr$SL_phase*24/(2*pi)))  
  Rhythmpcr=Rhythmpcr %>% add_column(LL_Phase_h = (Rhythmpcr$LL_phase*24/(2*pi))) 
  Rhythmpcr=Rhythmpcr %>% add_column(amp_difference = (Rhythmpcr$SL_amp-Rhythmpcr$LL_amp))
  Rhythmpcr=dplyr::select(Rhythmpcr,!c("SL_phase","LL_phase"))
  Rhythmpcr=mutate(Rhythmpcr, phase_difference=NA)
  for (i in 1:length(Rhythmpcr$SL_Phase_h)){
    Rhythmpcr$phase_difference[i] = hourCalculator(Rhythmpcr$SL_Phase_h[i],Rhythmpcr$LL_Phase_h[i])}
  return(Rhythmpcr)}
neuro_CR <- DiffRhythm_neuroSLvsLL("ST","mod_sel")

# plotting

qPCR_summary_diet_out %>% filter(Gene %in% neuro_genes)%>% 
  ggplot(aes(x=zt_on,y=LogFC,color=photoperiod,linetype=diet))+geom_point()+geom_linerange(aes(ymin=LogFC,ymax=LogFC+se))+geom_smooth(se=FALSE)+theme_classic()+facet_grid(Gene~diet,scale="free_y")+scale_color_manual(values=Photocolour)+scale_x_continuous(breaks=c(0,4,8,12,16,20,24))+labs(x="ZT")+theme(strip.text = element_text(face = "italic"))

neuro_box <- qPCR_tidy_out %>% filter(Gene %in% neuro_genes) %>% 
  ggplot(aes(x=diet,y=LogFC, color=photoperiod))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitterdodge(dodge.width = 0.75),alpha=0.5)+facet_wrap(~Gene, nrow = 1)+scale_color_manual(values=Photocolour)+theme_classic()+theme(strip.text = element_text(face = "italic"))

neuro_circ <- qPCR_tidy_out %>% filter(Gene %in% neuro_genes)%>% inner_join(neuro_CR,by = c("tissue_gene" = "id"))
neuro_circ$category <- factor(neuro_circ$category,levels=c("same","arrhy","loss","gain","change"))

neuro_circ_plot <- neuro_circ %>% ggplot(aes(x=ST,y=LogFC,color=photoperiod,linetype=category))+geom_point(aes(shape=diet),alpha=0.5)+cosinor_curve+theme_classic()+facet_wrap(~Gene,nrow = 1)+scale_color_manual(values=Photocolour)+scale_x_continuous(breaks=c(0,4,8,12,16,20,24))+labs(x="ST")+theme(strip.text = element_text(face = "italic"))

neuro_circ_plot/neuro_box

#2-way ANOVA stats

qPCR_tidy_out$tissue_gene <- factor(qPCR_tidy_out$tissue_gene)
df.aov <- qPCR_tidy_out %>%
  group_by(tissue_gene) %>%
  nest() %>%
  mutate(aov = map(data, ~aov(LogFC ~ photoperiod * diet, data = .x)))
df.aov$aov
pcr_res <- lapply(df.aov$aov,anova)
pcr_res <- t(data.frame(map(pcr_res, ~ select(., 5))))
rownames(pcr_res) <- df.aov$tissue_gene
pcr_res <- as.data.frame(pcr_res)

top_photo_pcr <- filter(pcr_res,photoperiod < 0.05)
qPCR_tidy_out %>% filter(tissue_gene %in% rownames(top_photo_pcr)) %>% 
  ggplot(aes(x=diet,y=LogFC, color=photoperiod))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitterdodge(dodge.width = 0.75),alpha=0.5)+facet_wrap(~Gene, scale="free_y",nrow = 4)+scale_color_manual(values=Photocolour)+theme_classic()
