# The specific lines of code that produces figures are labelled such as #Fig2B or #SFig1 or #DataS4.

library(data.table)
library(tidyverse)
library(MAPX)

dir.create("results")

customPlot <- list(
  theme_bw(base_size = 12),
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin=ggplot2::margin(5,5,5,5, "pt")
  )
)

#legend getting
get_legend<-function(myggplot) { # to extract legend from a single ggplot
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

Pal25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

reps <- c("0","1","2")

timepoints <- paste0("tp",c(4,10,16,22,28,34,40))
dir.create("results/normalization_plots")
datafiles <- list.files(path="/data/datafiles", pattern="Pf_Proteins.txt", full.names=TRUE)

all.features <- all.predictors <- all.fitdata <- all.scaledata <- all.rawdata <- vector(mode="list",length=7) %>% setNames(timepoints) %>%
  lapply(as.list)

for(df in datafiles) {
    
  # 1. Loading the raw data.
  
  rep <- str_extract(df,"(?<=IDC)[[:digit:]]+")
  tp <- paste0("tp",str_extract(df,"[[:digit:]]+(?=hpi)"))
  MCfile <- datafiles[grepl(paste0("IDC",rep,"_",str_remove(tp,"tp"),"hpi"),datafiles)]
  
  MCdata <- fread(MCfile) %>%
    setNames(make.names(names(.))) %>%
    dplyr::select(Accession,X..AAs,starts_with("Abundances..Grouped...")) %>%
    mutate(Accession=str_remove(Accession,".[[:digit:]]-p[[:digit:]]")) %>%
    left_join(MAPX::plasmoDB_data%>%dplyr::select(c(Accession,Product.Description))) %>%
    setNames(c("protein","length",paste0("Ab",1:10),"description")) %>%
    dplyr::select(protein,description,length,starts_with("Ab"))
  all.rawdata[[tp]][[rep]] <- MCdata

  # Clean up of raw data
  MCdata.clean <- X.clean.meltcurves(MCdata)
  
  if(nrow(MCdata.clean$data)==0) {
    next
  }
  # 3a. Extraction of raw data features
  features.raw <- X.extract.raw.features(MCdata.clean$data, abundances=c(1,2), remains=10)
  
  # 3b. Scaling of melting curves
  # The cleaned data are scaled from 0-1 and normalized into a sigmoid curve.
  MCdata.scaled <- X.scale.meltcurves(MCdata.clean$data)
  # The scaled data distribution before and after the normalization is plotted as a ggplot that
  # can be saved.
  ggsave(paste0("results/normalization_plots/MCdata_normalization_",tp,"_rep",rep,".png"),MCdata.scaled$plot)
  
  all.scaledata[[tp]][[rep]] <- MCdata.scaled$data

  # 4. Fitting of melting curves
  # For this function, the input data frame must contain the column 'protein' and one column for
  # each temperature starting with 'T' followed by a number.
  MCdata.fitted <- X.fit.meltcurves(MCdata.scaled$data)
  
  # The output gives a list of four elements.
  # The first element contains fitted values. Some of these can be used as features for complex prediction.
  
  # The fitdata will be important in later steps, let's save it in a list:
  all.fitdata[[tp]][[rep]] <- MCdata.fitted$data
  
  # The fitted melting curves can be visualized easily.
  #MCdata.curves <- X.plot.meltcurves(MCdata.scaled$data, MCdata.fitted$fits, pdf.name=paste0("MCcurves_",tp,"_rep",rep))
  
  # 5. Join raw features and fit features
  replace.inf <- function(x) {
    ifelse(is.infinite(x),
           ifelse(x < 0, min(x[is.finite(x)]-1, na.rm = TRUE), max(x[is.finite(x)]+1, na.rm = TRUE)),
           x
    )
  }
  
  features <- MCdata.fitted$data %>%
    # rownames_to_column("protein") %>%
    full_join(features.raw) %>%
    mutate(logABL=log2(ABL)) %>%
    mutate(logAUC=log10(AUC)) %>%
    mutate(logLoss=log10(Loss)) %>%
    mutate(Penalty=ifelse(is.na(Penalty),2,Penalty)) %>%
    mutate(Penalty_trans=log(Penalty)) %>%
    mutate(Penalty_trans=ifelse(is.infinite(Penalty_trans)&Penalty_trans<0,min(Penalty_trans[which(is.finite(Penalty_trans))]),Penalty_trans)) %>%
    mutate(across(where(is.numeric), ~ replace.inf(.))) %>%
    select(protein,Ti,logAUC,logABL,logLoss,Penalty_trans)
  
  all.features[[tp]][[rep]] <- features
  # The features define each single proteins. Relating the proteins to each other for each protein
  # will the calculation of the predictors.
  
  # 6. Calculation of the predictors
  # In this function, 'data' is a data frame that has to contain a column protein and one column
  # for each element in the vector given in 'features'. 'funs' are functions that will be applied
  # to calculate for each feature and 'prefixes' will be used to names the reuslting columns.
  # For example, in the code below, the first feature is 'Ti' and so the function applied to
  # 'Ti' between each pair of proteins will be absdif. Use ?X.calculate.predictors to see what
  # functions can be applied. The resulting column will be 'dTi' as the first prefix given
  # is 'd'.
  
  feature.predictors <- X.calculate.predictors(data=features,
                                               features=c("Ti","logAUC","logABL","logLoss","Penalty_trans"),
                                               funs=c(rep("absdif",4),"Xsum"),
                                               prefixes=c(rep("d",4),"sum")
  )
  
  # The previous function calculates 5 predictors from the features. An additional predictor is just a correlation
  # coefficient between data points of two proteins.
  
  cordata <- X.calculate.cors(MCdata.scaled$data)
  
  # The data frames feature.predictors and cordata should contain the same number of rows and rbind should be sufficient
  # to connect them. MAPX offers a generalized function to connect two data frames that contain two columns that are
  # interchangeable: MAPX::cross_join.
  
  all.predictors[[tp]][[rep]] <- MAPX::cross_join(cordata,feature.predictors,vars=paste0("protein",c(1,2)),mode="full")
    
}

saveRDS(all.predictors,"results/all.predictors.RDS")
saveRDS(all.features,"results/all.features.RDS")
saveRDS(all.fitdata,"results/all.fitdata.RDS")
saveRDS(all.scaledata,"results/all.scaledata.RDS")
saveRDS(all.rawdata,"results/all.rawdata.RDS")

## Assembly of the gold standard
#At this point, the data has been reduced into predictors that can be fed into the machine learning model. To train the model, gold standard dataset of interacting and non-interacting proteins has to be available.

### Assemble positive gold standard
#The positive gold standard is a pairwise list of protein-protein interactions that are already known. This, besides the data quality, is the most important model input, because the algorithm will search for protein interactions with similar properties as those of protein-protein interactions included in the positive gold standard. A care should be taken when putting this together: presumably, MAP-X only detect stable and stechiometrically defined interactions.Most databases list gold standard interactions with one column for complex name and one column for protein ID.The following command can take in this kind of data and convert it into a pair-wise interaction table.

GSpos <- X.assemble.posstandard(MAPX::gold.standard, sep=";")

### Assemble negative gold standard
#The most comprehensive approach to assemble the gold standard of non-interacting proteins is using GO terms: the proteins that do not interact cannot share GO terms. A care should be taken when predicting interactome for organisms without well annotated proteome - in such cases, using GO terms from a similar organism based on orthology would be a better strategy. The function expects a data frame with the first column containing protein IDs and the second column containing comma-separated GO terms.
#The following lines take a database from plasmoDB and prepares it for MAPX commands. Next, `X.assemble.negstandard` is used to assemble the standard of non-interacting proteins.

GOterms <-  MAPX::plasmoDB_data %>%
  dplyr::select(source_id,contains("GO", ignore.case=FALSE)) %>% # selects source_id column and all GO columns
  dplyr::select(contains("id")) %>% # selects all columns that contain id (gets rid of GO term descriptors)
  tidyr::unite("GO",contains("GO"),sep=";") %>% # connects all GO terms to semi-colon separate column "GO"
  #option 1: do not use any proteins that have N/A in at least one of the GO terms
  mutate(GO=ifelse(str_detect(GO,"N/A"),"",GO)) %>%
  #option 2: use everything. The algorithm will exclude those that have no GO term at all.
  # mutate(GO=ifelse(GO=="N/A;N/A;N/A","N/A",GO)) %>% # converts rows without GO to N/A
  # mutate(GO=str_remove_all(GO,"N/A;")) %>% # removes all N/As followed by semi-colon
  # mutate(GO=str_remove_all(GO,";N/A")) %>% # removes all N/As at the end of the row
  mutate(GO=ifelse(GO=="N/A","",GO)) %>%
  rename(protein=source_id) %>%
  mutate(protein=str_remove(protein,".[[:digit:]]$")) # removes '.1' at the end of the IDs

GSneg <- X.assemble.negstandard(data=GOterms, sep=";") %>%
  mutate(complex=0)

## Model training
#The positive and negative standards can be joined together

GS <- rbind(GSpos,GSneg)
fwrite(GS,"results/GS.csv")


all.predictors <- readRDS("results/all.predictors.RDS")
all.features <- readRDS("results/all.features.RDS")
all.fitdata <- readRDS("results/all.fitdata.RDS")
GS <- data.table::fread("results/GS.csv")


# Let's build a global protein-protein interaction network for tp28.
comp.models <- list()
mfolder <- "results/models"
dir.create(mfolder)

for(tp in timepoints) {
  for(rep in reps) {
    # 8. Model tuning
    # Machine learning models can be calculated with different model parameters. The best choice of the parameters
    # depends on the nature of the dataset. In this step, the best parameters for particular model can be chosen.
    # MAP-X comes with commands for tree-type models. The best-working model in my hands was random forest,
    # but any type of model could theoretically be integrated into the pipeline.
    # First, extract from gold standard (GS) those pairs of proteins that you have data from:
    
    ## IN THIS PERTICULAR CASE< THE MODELS WILL BE TRAINED ON ALL DATA ALL Pf-Hs PAIRS WILL BE EXCLUDED FROM PREDICTIONS (EXCEPT EXPORTED Pf PROTEINS)
    GS_specific <- X.assemble.traindata(all.predictors[[tp]][[rep]],GS)
    
    # This specific standard will be used to tune the parameters. 
    
    CV <- X.tune.tree(GS_specific,labels.col="complex",tree.type="RF",mTry=c(1,2,3,4,5),nodesize=c(5,10,15,20), downsample=3)
    
    # Let's look at the resulting data frame and choose the best model based on the chosen metric (in our case area under the precision-recall curve).
    # We can take the row with the best metric and extract the parameters as a named list - this named list is an argument for final model training.
    best.train.pars <- CV$data %>% arrange(desc(metric.mean)) %>% dplyr::select(!starts_with("metric")) %>% slice_head(n=3)
    # 9. X.crosstrain.tree is used for final model training. The function goes through a specified amount of train cycles. In each train cycle,
    # the data training data is randomly split to train.split number of subdata and one model is trained for each split. This results
    # in train.cycle x train.split number of models.
    
    for(btp in 1:nrow(best.train.pars)) {
      best.train.pars.row <- best.train.pars %>% slice(btp) %>% as.vector()
      cross.model <- X.crosstrain.tree(data=all.predictors[[tp]][[rep]],standard.set=GS_specific,train.pars=best.train.pars.row,
                                       evaluate=TRUE,plot=TRUE,train.split=4, train.cycles=3,
                                       save.dir=paste0(mfolder,"/crossmodels_",tp,"_rep",rep,"_btp",btp))
      # 10. The submodels can be averaged into a composite model using X.average.crossmodels
      comp.models <- X.average.crossmodels(data=cross.model$data, eval.metric="prc", evaluate=TRUE, standard.set=GS_specific)
      data.table::fwrite(comp.models$data, paste0(mfolder,"/crossmodels_",tp,"_rep",rep,"_btp",btp,"/averaged_predictions.csv"))
      ggsave(paste0(mfolder,"/crossmodels_",tp,"_rep",rep,"_btp",btp,"/averaged_predictions.png"), comp.models$plot, units="cm", height=8,width=8)
      rm(comp.models, cross.model)
      gc()
    }
  }
}

# For each timepoint and replicate, there are three different composite models. Precision-recall curves should be checked
# which models will be used in the next step for model averaging.
# An example file (model_choice.csv") is included in the github repository.

### Averaging predictions across replicates


nfolder <- "models/networks1"
dir.create(paste0("results/",nfolder))
# mchoice is a table with selected models from previous step
mchoice <- fread("/data/model_choice.csv") %>%
  setNames(c("rep",timepoints))

design=list('setNodeShapeDefault("ELLIPSE")',
            'setNodeColorDefault("#5A5A5A")',
            'lockNodeDimensions(TRUE)',
            'setNodeLabelOpacityDefault(1)',
            'setEdgeColorMapping(table.column="complex",list("1","0"), list("#00FF00","#FF0000"),mapping.type="d")',
            'setEdgeLineWidthMapping(mapping.type="c",table.column="score", widths=c(0.2,10))'
            
)

all.predictions <- list()
for(tp in timepoints) {
  rep.predictions <- list()
  for(rp in reps) {
    RFno <- mchoice %>% filter(rep==as.numeric(rp)) %>% pull(tp)
    if(is.na(RFno)) {
      next
    } else {
      rep.predictions[[rp]] <- fread(paste0(mfolder,"/crossmodels_",tp,"_rep",rp,"_btp",RFno,"/averaged_predictions.csv"))
    }
  }
  #averaging of predictions across replicates
  all.predictions[[tp]] <- X.average.reps(data=rep.predictions,
                                          evaluate=TRUE,standard.set=GS)
  
  ggsave(paste0("results/",nfolder,"/averaged.predictions_AUPRC_tp",tp,".png"), all.predictions[[tp]]$plot)
  fwrite((all.predictions[[tp]]$data),  paste0("results/",nfolder,"/averaged.predictions_tp",tp,".csv"))
  print(paste0("Done ", tp,"."))
}

for(tp in timepoints) {
  all.tp.predictions <- all.predictions[[tp]]
  
  statsplots <- list()
  statsdata <- list()
  final_networks <- list()
  gc()
  precision_scan <- seq(0.3,0.9,0.02)
  for(ps in precision_scan) {
    
    # lets find the closest score value to the current precision value
    cutoff=all.tp.predictions$eval.data %>%
      slice_min(abs(Precision-ps)) %>%
      pull(Probability) %>% unique() %>% .[1]
    # based on this value, the data will be cut
    cut_data <- all.tp.predictions$data %>% filter(score>=cutoff)
    
    # 12. The cut_data is prepared to be used for global network building. Let's first build the initial network:
    network.built <- X.build.complexes(data=cut_data, algo="SP",scores.col="score",init.stats=TRUE,final.stats=TRUE,
                                       standard.set=GS,labels.col="complex",  SP.finalpreds="original")
    
    # 13. The initial network is then refined.
    network.refined <- X.refine.complexes(data=network.built$data, algo="SCBP",final.stats=TRUE,
                                          standard.set=GS)
    
    # 14. Finally, the initial network is post-processed.
    network.post.split <- X.postprocess(data=network.refined$data, mode="split", scores.col="score", final.stats=TRUE,
                                        standard.set=GS,labels="complex", weighted=FALSE)
    
    network.post.trim <- X.postprocess(data=network.post.split$data, mode="trim", scores.col="score", final.stats=TRUE,
                                       standard.set=GS,labels="complex", weighted=TRUE)
    
    final_networks[[paste0("ps",ps)]] <- network.post.trim$data
    
    # Let's save the resulting plots and the network statistics into lists.
    
    statsplots[[paste0("ps",ps)]] <-
      network.built$stats_initial$plot + ggtitle("Averaged") +
      network.built$stats_final$plot + ggtitle("Built") +
      network.refined$stats_final$plot + ggtitle("Refined") +
      network.post.split$stats_final$plot + ggtitle("Split") +
      network.post.trim$stats_final$plot + ggtitle("Trimmed") +
      patchwork::plot_layout(ncol=5,nrow=1)
    
    statsdata[[paste0("ps",ps)]] <- network.post.trim$stats_final$data %>% bind_cols(network.post.trim$stats_final$network.stats) %>% mutate(PS=ps) %>%
      mutate()
    print(paste0("DONE WITH THE CUTOFF ", ps,"!"))
  }
  saveRDS(final_networks,paste0("results/",nfolder,"/",tp,"_final_networks.RDS"))
  
  # The following lines will calculate the network score from the statsdata. The network with cutoff with highest network score
  # is chosen as the final complexome map.
  
  choicedata <- statsdata %>%
    purrr::reduce(bind_rows) %>%
    mutate(network.score=stringency + fragmentation + purity + retrieval + recall/2)
  
  fwrite(choicedata, paste0("results/",nfolder,"/tp",tp,"_cutoff.choice.csv"))
}


# Plotting the complexome maps by Cytoscape
#Choosing the maps with the best network score, the global complexome map can be plotted in Cytoscape using the following function:
design=list('setNodeShapeDefault("ELLIPSE")',
            'setNodeColorDefault("#5A5A5A")',
            'lockNodeDimensions(TRUE)',
            'setNodeLabelOpacityDefault(1)',
            'setEdgeColorMapping(table.column="complex",list("1","0"), list("#00FF00","#FF0000"),mapping.type="d")',
            'setEdgeLineWidthMapping(mapping.type="c",table.column="score", widths=c(1,20))'
)

#final data here refers to the final_networks network with a selected cutoff.
# Fig1H #Fig2a #Fig2b #Fig2c

# This is commented out because it requires GUI with installed cytoscape to run
# final_networks <- list.files("/data/final_choices",full.names=TRUE)
# for(fn in final_networks) {
#   plotname = str_remove(fn,".csv")
#   networkfile <- fread(fn)
#   X.plot.network(networkfile,scores.col="pred",standard.set=GS,labels.col="complex",annotation=MAPX::plasmoDB_data%>%dplyr::rename(protein=Accession), 
#                  design.params=design, network.name="MAPX_network")
# }



# For local subnetworks, calibrated scores are typically used that approximate the probabilities
all.predictions.cal <- list()
for(tp in timepoints) {
  all.predictions.cal[[tp]] <- all.predictions[[tp]]$data %>%
    X.calibrate.scores(GS)
}

# Fig1D
all.predictions.cal <- all.predictions.cal %>%
  lapply(function(x) x$data) %>% 
  rbindlist(idcol="timepoint") %>%
  relocate(timepoint, .after=last_col())
fwrite(all.predictions.cal,"results/all.predictions.cal.csv")

all.predictions <- all.predictions %>%
  lapply(function(x) x$data) %>% 
  rbindlist(idcol="timepoint") %>%
  relocate(timepoint, .after=last_col())
fwrite(all.predictions,"results/all.predictions.csv")

# This is commented out because it requires GUI with installed cytoscape to run
#Calculation of local subnetwork of a protein of interest
#Fig3B #Fig3C #Fig3D #Fig3E #Fig3F #Fig3G #Fig2I #Fig3K (replace 'PF3D7_0412200' with the protein of interest)
# local.subnetwork <- X.local.network(all.predictions.cal%>%dplyr::select(!score),"PF3D7_0412200", plot="cytoscape",
#                                    plot.annotation=plasmoDB_data%>%rename(protein=Accession),
#                                    plot.cytoscape.path="C://Program Files/Cytoscape_v3.9.1/Cytoscape.exe")
# The command can be run alternatively with plot="visnetwork" which will export visnetwork graph
local.subnetwork <- X.local.network(all.predictions.cal%>%dplyr::select(!score),"PF3D7_1149100", 
                                    scores.col="probability", condition.col="timepoint", plot="visnetwork",
                                   plot.annotation=plasmoDB_data%>%rename(protein=Accession),
                                   plot.path="results/models/networks1")


# INTERACTION DYNAMICS (related to #Fig4)

all.features.data <- data.frame()
for(tp in timepoints) {
  for(rep in reps) {
    all.features.data <- bind_rows(all.features.data,
                                   all.features[[tp]][[rep]] %>% mutate(replicate=as.character(rep)) %>% mutate(condition=tp)
    )
  }
}
all.features.data <- all.features.data %>%
  filter(replicate!="0")

all.fitdata.frame <- data.frame()
for(tp in timepoints) {
  for(rep in reps) {
    all.fitdata.frame <- bind_rows(all.fitdata.frame,
                                   all.fitdata[[tp]][[rep]] %>% mutate(replicate=as.character(rep)) %>% mutate(condition=tp)
    )
  }
}
all.fitdata.frame <- all.fitdata.frame %>%
  filter(replicate!="0")

quality.data <- all.fitdata.frame %>%
  dplyr::select(protein,condition,replicate,R2) %>%
  rename(quality=R2) %>%
  filter(replicate!="0")

focus.complexes <- c('ribosome_full','proteasome','prefoldin','T_complex','proteasome_reg','E2_enzyme',
                     'microtubules','EF1','IF3','IF4F','CLAG_core','RNApolII_core','apicoplast_ribosome',
                     'HDP_complex','glyco_ENO_PK','histones_H2dimer','histones_H2Zdimer','VTA','Hsp_Hop',
                     'PROSC','IF2','spliceosome_snRNP_core','spliceosome_snRNP_U6_core','ETF1','histones_H3H4dimer','COPII_subcomplex1')

data.PCAs <- X.PCA.reduction(all.features.data,plot=TRUE,plot.complexes=MAPX::complexes[focus.complexes])

# The function X.AI will calculate the assembly factors.
# This takes a long time. Decreasing the number of trials (related to bootstrapping) will make it faster.
data.AIs <- X.AI(data.PCAs$data,MAPX::complexes[focus.complexes], quality.data, trials=50)

# Export of dynamics plots
dir.create("results/aveplots")
#Fig4C #Fig4D #Fig4E
for(cpx in names(data.AIs$plots)) {
  ggsave(paste0("results/aveplots/",cpx,".png"), data.AIs$plots[[cpx]])
}

# this table is available in the github repository
ctfu <- fread("/data/complexes_table_filledup.csv")

dir.create("results/Fig")

#Fig4B follows
plot <- data.AIs$data %>%
  filter(complex %in% focus.complexes) %>%
  mutate(condition=str_remove(condition,"tp")) %>%
  mutate(condition=factor(condition,levels=str_remove(timepoints,"tp"))) %>%
  dplyr::select(complex,condition,mean.AI) %>%
  distinct(.keep_all=TRUE) %>%
  group_by(complex) %>%
  mutate(average_z=scale(mean.AI)) %>%
  ungroup() %>%
  left_join(ctfu%>%dplyr::select(Complex,`Complex short name`)%>%setNames(c("complex","name"))%>%distinct(.keep_all=TRUE)) %>%
  ggplot(aes(x=condition,y=name,fill=average_z)) +
  geom_tile() +
  # scale_fill_gradient2(high="blue2",low="red", mid="black",na.value="gray80", limits=c(-3,3), midpoint=-0.25) +
  scale_fill_gradientn(colors=c("blue2","white","red","red"), values=scales::rescale(c(3,-0.25,-3,-8)), na.value="gray80") +
  scale_x_discrete(expand=c(0,0), name="Timepoint (hpi)") +
  scale_y_discrete(expand=c(0,0)) +
  customPlot +
  theme(legend.position="bottom",legend.title=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=4))

plotlegend <- get_legend(plot)
pdf("results/Fig/tileplot_z_averaged_centered_legend.pdf")
plot(plotlegend)
dev.off()
plot <- plot + theme(legend.position="none")
#Fig4B
ggsave(paste0("results/Fig/tileplot_z_averaged_centered.png"),plot,height=8,width=6,units="cm")

### Related to #Fig5 - moonlighting

my.rainbow <- c("red2","orange","yellow1","green3","cyan3","blue1","purple3")

table1_complexes <- MAPX::complexes[c("ribosome_full",  "IF2","IF3","EF1", "MSC_M","MSC_Q", "E2_enzyme", "apicoplast_ribosome", "SSU_UtpB", #protein turnover
                                      "proteasome", "proteasome_reg", # protein degradation
                                      "T_complex","R2TP","Hsp_Hop","prefoldin", # protein folding
                                      "MCM", "histones_H2Zdimer","histones_H2dimer","condensin1","condensin2","RuvBL_trimer",# chromatin regulation
                                      "RFC", # DNA replication
                                      "RNR","RNApolII", "NCBP","H_ACA", # RNA synthesis
                                      "spliceosome_snRNP_core","spliceosome_snRNP_U6_core", "spliceosome_U2_snRNP_S3Fa", # RNA splicing
                                      "qTRT", "KEOPS", "TRM", # RNA modification
                                      "RNH","exosome", #RNA degradation
                                      "SUMO_activating_enzyme", # post-translational modification
                                      "SRP","COPI","COPII","retromer","AP1", "PAM", # protein trafficking
                                      "NatA","glyco_ENO_PK","KDH","BCKDH","VTA", "PanK", # metabolism
                                      "microtubules", "basal_complex", #structural
                                      "PROSC", "HDP_complex", "CLAG_core","RON","RAP")] #plasmodium specific 

dir.create("results/moonlighters")

mapx_mls <- X.predict.moonlighters(data=data.PCAs$data,complexes=table1_complexes,min.reps=2,min.conditions=2,weights=data.PCAs$var.explained,
                                   kmeans.maxiter=100,kmeans.nstart=100, plot=TRUE)
#Fig5B
for(cpx in names(mapx_mls$PCA.plots)) {
  for(tp in str_remove(timepoints,"tp")) {
    for(rep in paste0("rep",reps)) {
      ggsave(paste0("results/moonlighters/",cpx,"_tp",tp,"_",rep,".png"),mapx_mls$PCA.plots[[cpx]][[tp]][[rep]])
    }
  }
}

S3_table <- mapx_mls$moonlighters %>%
  dplyr::select(complex,protein,condition,w.sig) %>%
  group_by(complex,protein) %>%
  summarise(Condition=paste(condition,collapse=", "), Significance=paste(w.sig,collapse=", ")) %>%
  ungroup() %>%
  setNames(c("Complex","Protein","Timepoints (hpi)","Significance")) %>%
  left_join(ctfu%>%dplyr::select(Complex,`Complex short name`)%>%distinct(.keep_all=TRUE)%>%setNames(c("Complex","Complex name")), by="Complex") %>%
  left_join(MAPX::plasmoDB_data%>%dplyr::select(Accession,Product.Description) %>%setNames(c("Protein","Protein name"))) %>%
  dplyr::select(starts_with("Complex"),starts_with("Protein"),starts_with("Timepoint"),Significance)
#DataS3
fwrite(S3_table,"results/moonlighters/DataS3_moonlighters.csv")

moonlighting_proteins_plot <- mapx_mls$moonlighters %>%
  # left_join(mapx_mls$data%>%dplyr::select(complex,protein,condition,),by=c("complex","protein","condition")) %>%
  mutate(condition=factor(condition,levels=gtools::mixedsort(unique(condition)))) %>%
  left_join(ctfu%>%dplyr::select(Complex,`Complex short name`)%>%setNames(c("complex","cname"))) %>%
  # mutate(cname=factor(cname,levels=intersect(unique(ctfu$`Complex short name`),unique(cname)))) %>%
  # mutate(significance=ifelse(significance<=0.05,significance,NA)) %>%
  ggplot(aes(y=protein,x=condition)) +
  ggforce::facet_col(facets = vars(cname), 
                     scales = "free_y", 
                     space = "free",
                     strip.position="top") +
  geom_line(aes(group=protein),size=0.1,color="gray50") +
  geom_point(aes(color=w.sig),size=1) +
  scale_color_gradientn(colours=c("brown","pink"),limits=c(0.000000005,0.05),na.value="gray50", transform="log10") +
  customPlot +
  theme_bw() +
  theme(axis.text.y=element_text(size=5,color="black"),
        axis.text.x=element_text(color="black"),
        panel.spacing = unit(0.1, "lines"),
        strip.text=element_text(size=7,hjust=0,margin=ggplot2::margin(1,0,0,0,"pt")),
        axis.title=element_blank(),
        strip.background=element_rect(color="transparent",fill="transparent"),
        # strip.text=element_blank(),
        legend.title=element_blank(),
        legend.position="bottom")
moonlighting_proteins_plot
#Fig5A
ggsave("results/moonlighters/moonlighting_proteins_all.png",moonlighting_proteins_plot,height=16,width=6,units="cm", dpi=600)

mapx_mls$moonlighters$complex %>% unique() %>% length()
mapx_mls$moonlighters$protein %>% unique() %>% length()
moonlighting_legend <- get_legend(moonlighting_proteins_plot)

pdf("results/moonlighters/moonlighting_legend_alt.pdf")
plot(moonlighting_legend)
dev.off()


#DataS4
PCAplots <- list()
for(cpx in names(table1_complexes)) {
  
  cpxn <- ctfu %>% filter(Complex==cpx) %>% pull(`Complex name`) %>% unique
  plotdata <- mapx_mls$data %>%
    dplyr::select(complex,protein,condition,cluster,replicate) %>%
    left_join(data.PCAs$data%>%dplyr::select(protein,condition,replicate,PC1,PC2),by=c("protein","condition","replicate")) %>%
    mutate(condition=factor(condition,levels=timepoints)) %>%
    mutate(replicate=as.character(as.numeric(replicate))) %>%
    filter(complex==cpx) %>%
    # filter(condition==tp) %>%
    # filter(replicate==rep) %>%
    na.omit()
  nosu <- plotdata %>% group_by(condition,replicate) %>% mutate(n=n_distinct(protein)) %>% ungroup() %>% pull(n) %>% max()
  labeldata <- plotdata %>% filter(cluster==2)
  if(nosu<3) {
    next
  }
  PCAplots[[cpx]] <- plotdata %>%
    ggplot(aes(x=PC1,y=PC2)) +
    facet_grid(condition~replicate) +
    ggforce::geom_mark_ellipse(aes(fill = factor(cluster)),size=0.1, expand=unit(5,"pt")) +
    ggrepel::geom_text_repel(data=labeldata, aes(x=PC1,y=PC2, label=protein), size=2, color="gray50") +
    geom_point(aes(color=as.character(cluster)), size=1.5) +
    scale_color_manual(values=c("black","red"),drop=FALSE) +
    scale_fill_manual(values=c("gray90","pink")) +
    scale_x_continuous(name="PC1", limits=c(-2.5,2.5)) +
    scale_y_continuous(name="PC2", limits=c(-2.5,2.5)) +
    customPlot +
    theme(legend.position="none",
          # axis.title=element_blank(),
          axis.text=element_text(size=7)) +
    ggtitle(cpxn)
}

#DataS4
pdf("results/moonlighters/moonlighting_all.pdf", width=8, height=15)
for (plot in PCAplots) {
  print(plot)
}
dev.off()



# Example meltcurves on one point
# alpha tubulin1 PF3D7_0903700
# beta tubulin PF3D7_1008700
# actin 1 PF3D7_1246200
# Fig1B

proteins <- list(Atub="PF3D7_0903700",Btub="PF3D7_1008700",act1="PF3D7_1246200") %>% stack() %>% setNames(c("protein","name"))
tofitdata <- all.scaledata$tp28$`1` %>%
  filter(protein %in% proteins$protein) %>%
  # filter(condition=="tp28" & replicate=="3") %>%
  dplyr::select(1:11)
  # column_to_rownames("protein") %>%
  # setnames(str_remove(names(.),"T"))

nowfits <- MAPX::X.fit.meltcurves(tofitdata)
nowfitted <- lapply(proteins$protein,
                    function(x) predict(nowfits$fits[[x]], newdata=data.frame(x=seq(37,73,length.out=100))) %>%
                      as.data.frame() %>% setNames("y") %>% add_column(x=seq(37,73,length.out=100))
) %>%
  setNames(proteins$protein) %>%
  bind_rows(.id="protein")

plot <- tofitdata %>%
  # rownames_to_column("protein") %>%
  pivot_longer(cols=!protein,names_to="x",values_to="y") %>%
  mutate(x=str_remove(x,"T")) %>%
  mutate(across(!protein, as.numeric)) %>%
  arrange(desc(protein)) %>%
  ggplot(aes(x=x,y=y,color=protein)) +
  geom_point(aes(color=protein)) +
  geom_line(data=nowfitted) +
  customPlot +
  scale_x_continuous(name="Temperature (°C)", limits=c(37,73), expand=c(0.05,0)) +
  scale_y_continuous(name="Fraction soluble",limits=c(0,1), expand=c(0.05,0)) +
  scale_color_manual(values=c("red","blue","green2")) +
  theme(
    legend.position="none",
    legend.title=element_blank()
  )
plot

#Fig1B
ggsave("results/Fig/3protein_example.png",plot,width=6,height=6,units="cm")


#Fig1C follows:
all.predictors <- all.predictors %>%
  lapply(function(x) x %>%
           lapply(function(y) y %>%
                    mutate(sumPenalty_trans_trans = log(sumPenalty_trans*(-1)+(-1)*min(sumPenalty_trans,na.rm=TRUE)+0.0000001))
           ))
boxplot.variables <- c("R2","dTi",'dlogAUC',"dlogABL","dlogLoss","sumPenalty_trans","sumPenalty_trans_trans")
boxplot.names <- c("R-Squared","ΔTi","Δlog(AUC)","Δlog(Abundance/AAs)","Δlog(Loss)","sum(Penalty)","log(Σ Penalty)")
prot3_boxplots_data <- all.predictors$tp28$`2` %>%
  # filter(if_all(starts_with("protein"), ~ . %in% proteins$protein)) %>%
  mutate(complex=0) %>%
  mutate(complex=ifelse(protein1=="PF3D7_0903700"&protein2=="PF3D7_1008700",1,
                        ifelse(protein1=="PF3D7_0903700"&protein2=="PF3D7_1246200",0,complex))) %>%
  na.omit() %>%
  dplyr::select(!c(sumPenalty_trans)) %>%
  pivot_longer(cols=!c(protein1,protein2,complex), names_to="Variable",values_to="Value") %>%
  mutate(int=ifelse(protein1=="PF3D7_0903700" & protein2 %in% proteins$protein,TRUE,FALSE))
prot3_boxplots_plot <- prot3_boxplots_data %>%
  ggplot(aes(x=as.character(complex),y=Value)) +
  geom_boxplot(outliers=FALSE) +
  geom_point(data= ~ filter(.x, int==TRUE), aes(color=as.character(complex)),size=3) +
  scale_color_manual(values=c("green2","blue")) +
  facet_wrap(~Variable,nrow=2, scales="free_y") +
  customPlot +
  theme(strip.text=element_blank(),
        panel.spacing.x=unit(0.8,"cm"),
        axis.title=element_blank(),
        axis.text.x=element_text(size=16, color="black"),
        axis.text.y=element_text(color="black"),
        legend.position="none") +
  scale_x_discrete(breaks=c(0,1),labels=c("-","+"))
#Fig1C
ggsave("results/Fig/3protein_example_predictors.png",units="cm",height=7,width=13)



#Fig1E follows:
# prediction statistics

all.predictions <- all.predictions %>%
  dplyr::rename(condition=timepoint) %>%
  dplyr::select(starts_with("protein"),condition,score)
dir.create("results/AUCs")
eval.metric.table <- data.frame()
prcplots <- list()
rocplots <- list()
for(tp in timepoints) {
  # load data
  predictions <- all.predictions %>% filter(condition==tp)
  
  toevaluate <- cross_join(predictions, GS, vars=paste0("protein",1:2), mode="inner")
  averaged_predictions.roc <- X.evaluate(toevaluate,scores.col="score",labels.col="complex", eval.metric="roc",plot=TRUE,black=TRUE)
  roc.plot <- averaged_predictions.roc %>%
    .$plot +
    theme(plot.title=element_blank(),
          axis.title=element_blank(),
          axis.text.x=element_text(size=20,color="black"),
          axis.text.y=element_blank(),
          legend.position="none") +
    scale_x_reverse(labels=c("0","0.25","0.5","0.75","1"), expand=c(0.001,0.001)) +
    scale_y_continuous(breaks=c(0.25,0.5,0.75),expand=c(0.001,0.001))
  
  averaged_predictions.prc <- X.evaluate(toevaluate,scores.col="score",labels.col="complex", eval.metric="prc",plot=TRUE,black=TRUE)
  prc.plot <- averaged_predictions.prc %>%
    .$plot +
    theme(plot.title=element_blank(),
          axis.title=element_blank(),
          axis.text.x=element_text(size=20,color="black"),
          axis.text.y=element_blank(),
          legend.position="none") +
    scale_x_continuous(labels=c("0","0.25","0.5","0.75","1"), expand=c(0.001,0.001)) +
    scale_y_continuous(breaks=c(0.25,0.5,0.75),expand=c(0.001,0.001))
  curdata <- data.frame(tp=tp,roc=averaged_predictions.roc$eval.metric,prc=averaged_predictions.prc$eval.metric)
  eval.metric.table <- bind_rows(eval.metric.table,curdata)
  
  #Fig1E
  ggsave(paste0("results/AUCs/AUROC_",tp,"hpi.png"),roc.plot,height=9, width=9, units="cm")
  ggsave(paste0("results/AUCs/AUPRC_",tp,"hpi.png"),prc.plot,height=9, width=9, units="cm")
  prcplots[[hpi]] <- prc.plot
  rocplots[[hpi]] <- roc.plot
  
}
#Fig1E
fwrite(eval.metric.table, "results/AUCs/eval.metric.table.csv")


#Fig1G follows:
# final.networks

# the folder final_choices contains csv files with interactions in complexome maps with highest network score
datafiles <- list.files(path="/data/final_choices", pattern="_pairwise.csv", full.names=TRUE) %>% gtools::mixedsort()
update_geom_defaults("text", list(size = 9))

all.stats <- list()
for (df in datafiles) {
  tp <- str_extract(df,"(?<=tp)[[:digit:]]+")
  cutoff <- str_extract(df,"(?<=cutoff)0\\.[[:digit:]]+")
  
  all.stats[[tp]] <- fread(df) %>% MAPX::X.network.stats(standard.set=GS,labels="complex",complex.stats=TRUE)
  
}

nowplots <- lapply(all.stats,
                   function(x) x$plot +
                     scale_y_continuous(limits=c(0,1500), name="No. of interactions") +
                     theme(
                       axis.text.x=element_text(size=24,color="black"),
                       axis.text.y=element_blank(),
                       axis.title=element_blank()
                     ))
#Fig1G
for(tp in timepoints) {
  ggsave(paste0("results/Fig/all.stats_tp",tp,".png"),nowplots[[as.character(tp)]], units="cm",height=9, width=9)
}


#Fig4A below:
# explanatory for timepoint analysis

plot <- data.PCAs$data %>% 
  filter(sample=="tp34;_;2") %>%
  mutate(complex=ifelse(protein %in% MAPX::complexes$T_complex,"YES","NO")) %>%
  arrange(complex) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=complex,size=complex,alpha=complex)) +
  scale_size_manual(values=c(0.25,1.5)) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_color_manual(values=c("Gray80","red")) +
  customPlot +
  scale_x_continuous(limits=c(-2.5,2.5), name=paste0("PC1 (explains ",data.PCAs$var.explained[1],"% of variation)")) +
  scale_y_continuous(limits=c(-2.5,2.5), name=paste0("PC2 (explains ",data.PCAs$var.explained[2],"% of variation)")) +
  theme(legend.position="none")
ggsave("results/Fig/PCA_plot_Tcomplex_rep3_34hpi.png", plot, units="cm",height=8,width=8)

# explanatory for timepoint analysis

plot <- all.features.data %>%
  filter(condition=="tp34") %>%
  filter(replicate=="2") %>%
  mutate(complex=ifelse(protein %in% MAPX::complexes$T_complex,"YES","NO")) %>%
  arrange(complex) %>%
  dplyr::select(!c(condition,replicate)) %>%
  mutate(across(!c(protein,complex),scale)) %>%
  pivot_longer(cols=!c(protein,complex),names_to="predictor",values_to="value") %>%
  mutate(predictor=factor(predictor,levels=c("Ti","logAUC","logABL","logLoss","Penalty_trans"))) %>%
  ggplot(aes(x=predictor,y=value)) +
  geom_boxplot(aes(fill=complex,alpha=complex, color=complex)) +
  scale_fill_manual(values=c("gray80","red")) +
  scale_color_manual(values=c("gray80","red")) +
  scale_alpha_manual(values=c(0.5,1)) +
  customPlot +
  scale_y_continuous(name="Scaled predictor value") +
  theme(legend.position="none",
        axis.title.x=element_blank())

ggsave("results/Fig/features_plot_Tcomplex_rep3_34hpi.png", plot, units="cm",height=8,width=8)

#No. of proteins detected
#SFig1B inlet
plot_no.proteins <- all.rawdata %>%
  lapply(function(x) x %>% rbindlist(idcol="replicate")) %>%
  rbindlist(idcol="condition") %>%
  mutate(condition=factor(condition,levels=timepoints)) %>%
  group_by(replicate,condition) %>%
  summarise(No.proteins=n()) %>%
  ungroup() %>%
  ggplot(aes(x=condition,y=No.proteins)) +
  # geom_boxplot(outlier.shape=NA) + 
  geom_line(aes(group=replicate)) +
  geom_point(aes(shape=replicate,color=condition),size=3)  + 
  customPlot +
  scale_x_discrete(name="hpi") +
  scale_y_continuous(name="Proteins detected")  +
  scale_color_manual(values=my.rainbow) +
  theme(plot.title=element_text(size=18),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.position="none")
plot_no.proteins  
ggsave(paste0("results/Fig/plot_no.proteins.png"),height=9, width=9, units="cm")

#all melting curves vs selected complexes
#SFig1A
plotcomplexes <- c("ribosome_full","T_complex","IF3","proteasome","proteasome_reg")
for(pc in plotcomplexes) {
  plot <- all.scaledata$tp34$`2` %>%
    dplyr::select(protein, matches("T[[:digit:]]+")) %>%
    pivot_longer(cols=!protein,names_to="Temperature",values_to="Scaled_abundance") %>%
    mutate(Temperature=as.numeric(str_remove(Temperature,"T"))) %>%
    mutate(ribosome=ifelse(protein %in% MAPX::complexes[[pc]],TRUE,FALSE)) %>%
    arrange(ribosome) %>%
    ggplot(aes(x=Temperature,y=Scaled_abundance)) +
    geom_jitter(aes(color=ribosome,size=ribosome)) +
    scale_color_manual(values=c("gray","red")) +
    scale_size_manual(values=c(0.1,0.5)) +
    scale_x_continuous(name="Melting temperature (°C)") +
    scale_y_continuous(name="Scaled abundance") +
    customPlot +
    theme(legend.position="none",
          axis.text.y=element_text(size=6, color="black",margin=ggplot2::margin(r=1)),
          axis.text.x=element_text(size=6, color="black",vjust=3),
          axis.title.y=element_text(size=7, color="black",margin=ggplot2::margin(0,0,0,0),vjust=2),
          axis.title.x=element_text(size=7, color="black",margin=ggplot2::margin(0,0,0,0),vjust=2.5),
          panel.border = element_rect(color="black",fill=NA,size=1),
          axis.ticks.length=unit(.025, "cm"),
          axis.ticks=element_line(linewidth=0.25))
  
  ggsave(paste0("results/Fig/allcurves_jitterplot_",pc,"_tp34_rep3.png"), plot, units="cm", height=3.5,width=3.5)
}

#Fig1A
plotpcomplexes <- MAPX::complexes[plotcomplexes] %>% stack() %>% setNames(c("protein","complex"))
plot <- all.scaledata$tp34$`2` %>%
  dplyr::select(protein, matches("T[[:digit:]]+")) %>%
  pivot_longer(cols=!protein,names_to="Temperature",values_to="Scaled_abundance") %>%
  mutate(Temperature=as.numeric(str_remove(Temperature,"T"))) %>%
  left_join(plotpcomplexes) %>%
  mutate(size=ifelse(is.na(complex), "small","big")) %>%
  arrange(desc(size)) %>%
  ggplot(aes(x=Temperature,y=Scaled_abundance)) +
  geom_jitter(aes(color=complex,size=size)) +
  scale_color_manual(values=c("red","blue1","green2","cyan2","magenta"), na.value="gray80") +
  scale_size_manual(values=c(1.5,0.6)) +
  scale_x_continuous(name="Melting temperature (°C)") +
  scale_y_continuous(name="Scaled abundance") +
  customPlot +
  theme(
    legend.position="none",
    axis.text=element_text(size=16),
    axis.title=element_text(size=18)
  )
ggsave(paste0("results/Fig/allcurves_jitterplot_mixed_tp34_rep3.png"), plot, units="cm", height=12,width=12)

#Fig2a #Fig2b below:
dir.create("results/Fig/individual_complexes")
datafiles <- list.files(path="/data/final_choices", pattern="_pairwise.csv", full.names=TRUE) %>% gtools::mixedsort()
all.maps <- list()
for (df in datafiles) {
  tp <- str_extract(df,"tp[[:digit:]]+")
  cutoff <- str_extract(df,"(?<=cutoff)0\\.[[:digit:]]+")
  
  all.maps[[tp]] <- fread(df)
  
}

complexes <- MAPX::complexes
complexes$MutS_DNApol <- c("PF3D7_1427500","PF3D7_0505500","PF3D7_1017000")
complexes$exosome_ring <- c("PF3D7_1427800","PF3D7_1340100","PF3D7_1364500","PF3D7_0209200","PF3D7_0412200","PF3D7_0626200")
plotcomplexes <- list(ribosome_full=40,
                      T_complex=34,
                      proteasome_reg=34,
                      IF3=28,
                      IF2=40,
                      microtubules=28,
                      CKII=34,
                      CLAG_core=4,
                      RuvBL_trimer=34,
                      EF1=40,
                      qTRT=34,
                      Hsp_Hop=22,
                      glyco_ENO_PK=10,
                      MSC_core=34,
                      KDH=4,
                      RNR=34,
                      exosome=40,
                      spliceosome_snRNP_U6_core=4,
                      RON=4,
                      MutS_DNApol=40,
                      exosome_ring=34
) %>% 
  stack() %>% setNames(c("timepoint","complex")) %>%
  dplyr::select(complex,timepoint)
all.scaledata.frame <- all.scaledata %>%
  lapply(function(x) x %>% rbindlist(idcol="replicate")) %>%
  rbindlist(idcol="condition")
for(i in 1:nrow(plotcomplexes)) {
  
  nowcomplex <- plotcomplexes %>% as_tibble() %>% dplyr::slice(i) %>% pull(complex) %>% as.character()
  proteins <- complexes[[nowcomplex]]
  tp <- plotcomplexes %>% dplyr::slice(i) %>% pull(timepoint) %>% paste0("tp",.)
  
  nowtable <- all.maps[[tp]] %>%
    filter(if_any(starts_with("protein"), ~ . %in% proteins))
  complexlist <- MAPX::X.pairwise.to.complexes(data=nowtable,scores.col="pred")
  # gold_proteins <- gold_table %>% dplyr::select(starts_with("protein")) %>% unlist() %>% unname() %>% unique()
  
  for(cl in seq_along(complexlist)) {
    listproteins <- complexlist[[cl]]
    listmembers <- intersect(listproteins, complexes[[nowcomplex]])
    
    for(rep in reps) {
      
      # meltcurves
      tofitdata <- all.scaledata.frame %>%
        filter(protein %in% listproteins) %>%
        filter(replicate==rep) %>%
        filter(condition==paste0("tp",plotcomplexes$timepoint[i])) %>%
        dplyr::select(3:13) 
      
      if(nrow(tofitdata)<2) {
        next
      }
      
      nowfits <- MAPX::X.fit.meltcurves(tofitdata)
      nowfitted <- lapply(names(nowfits$fits),
                          function(x) predict(nowfits$fits[[x]], newdata=data.frame(x=seq(37,73,length.out=100))) %>%
                            as.data.frame() %>% setNames("y") %>% add_column(x=seq(37,73,length.out=100))
      ) %>%
        setNames(names(nowfits$fits)) %>%
        bind_rows(.id="protein") %>%
        mutate(member=ifelse(protein %in% listmembers,TRUE,FALSE)) %>%
        mutate(member=factor(member,levels=c(FALSE,TRUE)))
      
      plot <- tofitdata %>%
        # rownames_to_column("protein") %>%
        mutate(member=ifelse(protein %in% listmembers,TRUE,FALSE)) %>%
        mutate(member=factor(member,levels=c(FALSE,TRUE))) %>%
        pivot_longer(cols=!c(protein,member),names_to="x",values_to="y") %>%
        mutate(x=str_remove(x,"T")) %>%
        mutate(across(!c(protein,member), as.numeric)) %>%
        ggplot(aes(x=x,y=y,color=member)) +
        geom_line(data=nowfitted, aes(group=protein), linewidth=0.15) +
        geom_point(aes(fill=member),size=1.8,shape=21,color="black",stroke=0.5) +
        customPlot +
        scale_x_continuous(name="Temperature (°C)", limits=c(37,73), expand=c(0.05,0)) +
        scale_y_continuous(name="Fraction soluble",limits=c(-0.03,1.03), expand=c(0.05,0)) +
        scale_color_manual(values=c("pink","cyan2"), drop=FALSE) +
        scale_fill_manual(values=c("pink","cyan2"), drop=FALSE) +
        theme(
          legend.position="none",
          legend.title=element_blank(),
          axis.text.y=element_text(size=9,color="black"),
          axis.text.x=element_blank(),
          axis.title=element_blank()
        )
      #Fig2A #Fig2B
      ggsave(paste0("results/Fig/individual_complexes/",nowcomplex,"_cl",cl,"_",tp,"_rep",rep,"_meltcurves.png"),plot,width=4,height=4,units="cm")
    }
  }
}

#SFig1D
# Histogram of sizes of clusters at different timepoints
plot <- all.maps %>%
  lapply(function(x) x %>% X.pairwise.to.complexes(scores.col="pred")) %>%
  lapply(function(x) lapply(x, function(y) y%>%length())) %>% 
  lapply(function(x) x%>%unlist%>%sort(decreasing=TRUE)%>%table%>%as.data.frame%>%setNames(c("Size","Frequency"))) %>% 
  data.table::rbindlist(idcol="Timepoint") %>% 
  mutate(Size=as.numeric(as.character(Size))) %>%
  mutate(Timepoint=factor(Timepoint, levels=gtools::mixedsort(unique(Timepoint)))) %>%
  ggplot(aes(x=Size,y=Frequency,color=Timepoint)) + 
  geom_line(aes(group=Timepoint)) + geom_point() +
  scale_color_manual(values=my.rainbow, name="Timepoint (hpi)") + 
  scale_x_continuous(trans="log2", name="Cluster size") +
  scale_y_continuous(trans="log2", name="No. of clusters") +
  theme_bw() +
  customPlot
ggsave("results/Fig/histogram_size_of_clusters.png", plot, units="cm", height=9, width=12)



# Example meltcurves for prefoldin at 28 hpi
#SFig2B
proteins <- complexes$prefoldin %>% as.data.frame() %>% setNames("protein")

tofitdata <- all.scaledata.frame %>%
  filter(protein %in% proteins$protein) %>%
  filter(condition=="tp28" & replicate=="0") %>%
  dplyr::select(3:13)

nowfits <- MAPX::X.fit.meltcurves(tofitdata)
nowfitted <- lapply(proteins$protein,
                    function(x) predict(nowfits$fits[[x]], newdata=data.frame(x=seq(37,73,length.out=100))) %>%
                      as.data.frame() %>% setNames("y") %>% add_column(x=seq(37,73,length.out=100))
) %>%
  setNames(proteins$protein) %>%
  bind_rows(.id="protein")


plot <- tofitdata %>%
  pivot_longer(cols=!protein,names_to="x",values_to="y") %>%
  mutate(x=str_remove(x,"T")) %>%
  mutate(across(!protein, as.numeric)) %>%
  arrange(desc(protein)) %>%
  ggplot(aes(x=x,y=y,color=protein)) +
  geom_point(aes(color=protein)) +
  geom_line(data=nowfitted) +
  customPlot +
  scale_x_continuous(name="Temperature (°C)", limits=c(37,73), expand=c(0.05,0)) +
  scale_y_continuous(name="Fraction soluble",limits=c(0,1), expand=c(0.05,0)) +
  scale_color_manual(values=Pal25) +
  theme(
    legend.position="right",
    legend.title=element_blank()
  )
#SFig2B
ggsave("results/Fig/prefoldin_curves.png",plot,width=10,height=6,units="cm")


# BELOW GENERATION OF CHARTS FOR FIG3:

complexes <- list(HDP=MAPX::complexes$HDP_complex,
                  NCBP=MAPX::complexes$NCBP,
                  PTP=c("PF3D7_0301700","PF3D7_0202200"),
                  SEMP=c('PF3D7_0702400','PF3D7_1353100','PF3D7_1126200','PF3D7_1026800','PF3D7_1003500','PF3D7_0705700','PF3D7_1465900','PF3D7_0322900','PF3D7_1105400','PF3D7_1342000','PF3D7_0517000','PF3D7_1004000','PF3D7_1431700','PF3D7_1351400','PF3D7_1341200','PF3D7_0618300','PF3D7_1142500','PF3D7_1027800','PF3D7_1142600','PF3D7_1109900','PF3D7_0706400','PF3D7_0304400','PF3D7_1424100','PF3D7_1323100','PF3D7_0714000'),
                  Alba=c("PF3D7_0814200","PF3D7_1006200","PF3D7_1346300","PF3D7_1347500"),
                  proteasome_reg_ub=c("PF3D7_0205900","PF3D7_0305700","PF3D7_0312300","PF3D7_0413600","PF3D7_0527100",
                                      "PF3D7_0527200","PF3D7_1008400","PF3D7_1129200","PF3D7_1130400","PF3D7_1225800",
                                      "PF3D7_1311500","PF3D7_1338100","PF3D7_1402300"),
                  RMC=c("PF3D7_1352800","PF3D7_1107200","PF3D7_1226400"),
                  exp=c("PF3D7_0301700","PF3D7_1201200","PF3D7_1149100"))

#Detection data
detection_table <- all.rawdata %>%
  lapply(function(x) x %>% rbindlist(idcol="replicate")) %>%
  rbindlist(idcol="condition") %>%
  dplyr::select(!replicate) %>%
  distinct(.keep_all=TRUE) %>%
  mutate(timepoint=as.integer(str_remove(condition,"tp"))) %>%
  dplyr::rename(subunit=protein) %>%
  left_join(stack(complexes)%>%setNames(c("subunit","complex")), relationship="many-to-many") %>%
  dplyr::select(subunit,timepoint,complex) %>% 
  mutate(n_codetected=-1) %>% mutate(perc_detected=-1)

# Global maps
all_complexes <- list()
all_networks_numbered <- list()
for(i in names(all.maps)) {
  all_complexes[[i]] <- X.pairwise.to.complexes(all.maps[[i]])
  all_networks_numbered[[i]] <- X.complexes.to.pairwise(all_complexes[[i]]) %>%
    mutate(timepoint=i)
  cat("Done", i,"\n")
}
all_networks_numbered <- reduce(all_networks_numbered,bind_rows)
# all_networks <- reduce(all_networks,bind_rows)

table1 <- list()
# detection in final networks
for(complex in names(complexes)) {
  print(complex)
  
  for(subunit in complexes[[complex]]) {
    table1[[paste0(complex,";;",subunit)]] <- all_networks_numbered %>% 
      filter(protein1==subunit | protein2==subunit) %>%
      pivot_longer(cols=c(protein1,protein2),values_to="protein",names_to=NULL) %>%
      filter(protein!=subunit) %>%
      mutate(isincomplex=ifelse(protein %in% complexes[[complex]],1,0)) %>%
      group_by(timepoint) %>%
      summarize(n_codetected=sum(isincomplex)) %>%
      ungroup() %>%
      mutate(perc_detected=n_codetected/length(complexes[[complex]])) %>%
      mutate(complex=complex) %>%
      mutate(subunit=subunit)
  }
}
table1_red <- purrr::reduce(table1, bind_rows)

#local subnetworks
all_proteins <- detection_table$subunit %>% unique()
table1.1 <- list()

# detection in local networks
for(complex in names(complexes)) {
  print(complex)
  
  for(subunit in complexes[[complex]]) {
    
    pretable <- all.predictions.cal %>%
      dplyr::rename(condition=timepoint,calpred=probability) %>%
      dplyr::select(protein1,protein2,condition,calpred) %>%
      mutate(condition=paste0("tp",condition)) %>% 
      X.local.network(proteins=subunit,min.conditions=1,plot=NULL, cutoff=0.4, condition.col="condition")
    
    if(nrow(pretable)>0) {
      table1.1[[paste0(complex,";;",subunit)]] <- pretable %>%
        filter(protein1==subunit | protein2==subunit) %>%
        dplyr::rename(timepoint=condition) %>%
        pivot_longer(cols=c(protein1,protein2),values_to="protein",names_to=NULL) %>%
        filter(protein!=subunit) %>%
        distinct(.keep_all=TRUE) %>%
        mutate(isincomplex=ifelse(protein %in% complexes[[complex]],1,0)) %>%
        group_by(timepoint) %>%
        summarize(n_codetected=sum(isincomplex)) %>%
        ungroup() %>%
        mutate(perc_detected=n_codetected/length(complexes[[complex]])) %>%
        mutate(complex=complex) %>%
        mutate(subunit=subunit)
      
    }
    
  }
}
table1.1_red <- purrr::reduce(table1.1, bind_rows)

library(patchwork)
for(cpx in names(complexes)) { 
  tileplot_global <- table1_red %>%
    mutate(timepoint=as.integer(str_remove(timepoint,"tp"))) %>%
    bind_rows(detection_table) %>%
    distinct(timepoint,complex,subunit, .keep_all=TRUE) %>%
    mutate(timepoint=factor(timepoint,levels=str_remove(timepoints,"tp"))) %>%
    group_by(complex) %>%
    complete(subunit, timepoint,fill=list(perc_detected=NA)) %>%
    ungroup() %>%
    filter(complex==cpx) %>%
    ggplot(aes(x=timepoint,y=subunit,fill=perc_detected)) +
    geom_tile(color="gray90",size=0.5) +
    scale_x_discrete(drop=FALSE,expand=c(0,0)) +
    scale_fill_gradientn(colors=c("gray50","red","orange2","yellow2","green3"), values = scales::rescale(c(-1,-0.99,0,0.0002, 0.05, 0.2, 0.55, 1)), limits=c(-1,1),na.value="white") +
    coord_fixed() +
    customPlot +
    theme(
      panel.background = element_rect(fill = "transparent", colour = "white"),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(color="gray90",fill=NA,size=unit(0,"pt")),
      plot.margin = unit(c(0,0,0,0), "pt"),
      panel.spacing=unit(-10,"lines"),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(margin=ggplot2::margin(1,0,0,0,"pt"),size=7, color="black"),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position="none")    
  #tileplot_global
  
  
  tileplot_local <- table1.1_red %>%
    mutate(timepoint=as.integer(str_remove(timepoint,"tp"))) %>%
    bind_rows(detection_table) %>%
    distinct(timepoint,complex,subunit, .keep_all=TRUE) %>%
    mutate(timepoint=factor(timepoint,levels=str_remove(timepoints,"tp"))) %>%
    group_by(complex) %>%
    complete(subunit, timepoint,fill=list(perc_detected=NA)) %>%
    ungroup() %>%
    filter(complex==cpx) %>%
    ggplot(aes(x=timepoint,y=subunit,fill=perc_detected)) +
    geom_tile(color="gray90",size=0.5) +
    scale_x_discrete(drop=FALSE,expand=c(0,0)) +
    scale_fill_gradientn(colors=c("gray50","red","orange2","yellow2","green3"), values = scales::rescale(c(-1,-0.99,0,0.0002, 0.05, 0.2, 0.55, 1)), limits=c(-1,1),na.value="white") +
    coord_fixed() +
    customPlot +
    theme(
      panel.background = element_rect(fill = "transparent", colour = "gray90"),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(color="gray90",fill=NA,size=unit(0,"pt")),
      plot.margin = unit(c(0,0,0,0), "pt"),
      panel.spacing=unit(-10,"lines"),
      axis.text.x = element_text(margin=ggplot2::margin(1,0,0,0,"pt"),size=7, color="black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = "none")
  
  #tileplot_local
  
  #library(patchwork)
  tileplot <- tileplot_global + plot_spacer() + tileplot_local + plot_layout(widths=c(1,0.005,1))
  ggsave(paste0("results/Fig/tileplot_",cpx,".png"),width=9, units="cm", dpi=600)

}

# Sequence alignments
#SFig2C #SFig2D #SFig4

# not run in the codeocean environment (problem setting up the environment for package msa)
# 
# library(msa)
# fastafiles <- paste0("Fig1/PfAlba",c(1,2,3,4),".fasta")
# for(ff in fastafiles) {
#   prot=str_extract(ff,"Alba[[:digit:]]")
#   sequences <- readAAStringSet(ff)
#   alignment <- msa(sequences)
#   msaPrettyPrint(alignment, output="tex", showNames="none", paperWidth=7.5,
#                  showLogo="none", askForOverwrite=FALSE, verbose=FALSE,
#                  shadingMode="similar",shadingColors="reds",
#                  showConsensus="none")
#   texfile <- readLines("alignment.tex")
#   texfile <- sapply(texfile, function(x) 
#     str_replace(x,"SAMUEL~1.PAZ","samuel.pazicky")
#   ) %>% unname()
#   writeLines(texfile,"alignment.tex")
#   tools::texi2pdf("alignment.tex", clean=TRUE)
#   file.rename("alignment.pdf","alignment_Rrp46.pdf")
#   saveWidth <- getOption("width")
#   options(width=100)
#   sink("myAlignment.txt")
#   print(alignment, show="complete", halfNrow=-1)
#   sink()
#   options(width=saveWidth)
#   file.rename("myAlignment.txt","alignment_Rrp46.txt")
# }
# 
# # Sequence alignment of PF3D7_0412200 and ScRRP46
# 
# library(msa)
# fastafiles <- paste0("Fig1/Rrp46.fasta")
# for(ff in fastafiles) {
#   sequences <- readAAStringSet(ff)
#   PAL <- pairwiseAlignment(pattern = sequences[[1]], subject = sequences[[2]], substitutionMatrix="BLOSUM50")
#   wPAL <- c(
#     paste0(">",names(sequences[1])),
#     PAL@pattern %>% as.character(),
#     paste0(">",names(sequences[2])),
#     PAL@subject %>% as.character()
#   )
#   wPALname <- str_split_1(fastafiles,"/") %>% .[length(.)] %>% str_remove(".fasta")
#   writeLines(wPAL,paste0("pairwiseAlignment_",wPALname,".txt"))
#   alignment <- msa(sequences)
#   msaPrettyPrint(alignment, output="tex", showNames="none", paperWidth=7.5,
#                  showLogo="none", askForOverwrite=FALSE, verbose=FALSE,
#                  shadingMode="similar",shadingColors="reds",
#                  showConsensus="none")
#   texfile <- readLines("alignment.tex")
#   texfile <- sapply(texfile, function(x) 
#     str_replace(x,"SAMUEL~1.PAZ","samuel.pazicky")
#   ) %>% unname()
#   writeLines(texfile,"alignment.tex")
#   tools::texi2pdf("alignment.tex", clean=TRUE)
#   file.rename("alignment.pdf","alignment_Rrp46.pdf")
#   saveWidth <- getOption("width")
#   options(width=100)
#   sink("myAlignment.txt")
#   print(alignment, show="complete", halfNrow=-1)
#   sink()
#   options(width=saveWidth)
#   file.rename("myAlignment.txt","alignment_Rrp46.txt")
#   
# }
# 
# # Sequence alignment of PF3D7_0626200 and ScRRP43
# 
# library(msa)
# fastafiles <- paste0("Fig/Rrp43.fasta")
# for(ff in fastafiles) {
#   sequences <- readAAStringSet(ff)
#   PAL <- pairwiseAlignment(pattern = sequences[[1]], subject = sequences[[2]], substitutionMatrix="BLOSUM50")
#   wPAL <- c(
#     paste0(">",names(sequences[1])),
#     PAL@pattern %>% as.character(),
#     paste0(">",names(sequences[2])),
#     PAL@subject %>% as.character()
#   )
#   wPALname <- str_split_1(fastafiles,"/") %>% .[length(.)] %>% str_remove(".fasta")
#   writeLines(wPAL,paste0("pairwiseAlignment_",wPALname,".txt"))
#   
#   
#   alignment <- msa(sequences)
#   msaPrettyPrint(alignment, output="tex", showNames="none", paperWidth=7.5,
#                  showLogo="none", askForOverwrite=FALSE, verbose=FALSE,
#                  shadingMode="similar",shadingColors="reds",
#                  showConsensus="none")
#   texfile <- readLines("alignment.tex")
#   texfile <- sapply(texfile, function(x) 
#     str_replace(x,"SAMUEL~1.PAZ","samuel.pazicky")
#   ) %>% unname()
#   writeLines(texfile,"alignment.tex")
#   tools::texi2pdf("alignment.tex", clean=TRUE)
#   file.rename("alignment.pdf","alignment_Rrp43.pdf")
#   
#   saveWidth <- getOption("width")
#   options(width=100)
#   sink("myAlignment.txt")
#   print(alignment, show="complete", halfNrow=-1)
#   sink()
#   options(width=saveWidth)
#   file.rename("myAlignment.txt","alignment_Rrp43.txt")
#   
# }
# 
# detach("package:msa", unload=TRUE)
# detach("package:Biostrings", unload=TRUE)
# detach("package:GenomeInfoDb", unload=TRUE)
# detach("package:XVector", unload=TRUE)
# detach("package:IRanges", unload=TRUE)