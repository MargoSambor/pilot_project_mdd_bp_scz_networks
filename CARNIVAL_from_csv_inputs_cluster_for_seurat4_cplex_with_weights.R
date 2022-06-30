library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(ggrepel)

library(DESeq2)

library(vsn)

library(reshape)
library(gridExtra)

library(cowplot)
library(hexbin)

library(limma)
library("purrr")

library("GEOquery")
library("pasilla")
#library("IHW")

# We will need this so we can use the pipe: %>%
library(magrittr)
library(tidyverse)

library("biomaRt")
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
library(igraph)
library(clusterProfiler)

#for carnival
library(CARNIVAL)
library(OmnipathR)
library(visNetwork)

## We also load the support functions for carnival
source_path="/nfs/research/petsalaki/users/sambor/pilot_project_mdd_bp_scz_networks-main/"
source(paste0(source_path, "assignPROGENyScores.r"))
source(paste0(source_path, "generateTFList.r"))
source(paste0(source_path, "carnival_visNetwork.r"))
#support functions for dorothea:
source(paste0(source_path, "support_functions.R"))

#tutorial on dorothea here:
#https://github.com/saezlab/transcriptutorial/blob/master/scripts/04_TranscriptionFactor_activity_with_Dorothea.md
#tutorial on preparing the normalized counts for the dorothea input here:
#https://github.com/saezlab/transcriptutorial/blob/master/scripts/01_normalisation.md
#tutorial on differential analysis for dorothea input here:
#https://github.com/saezlab/transcriptutorial/blob/master/scripts/02_differential_analysis.md
#CARNIVAL home:
#https://saezlab.github.io/CARNIVAL/
#CARNIVAL github tutorial:
#https://github.com/saezlab/transcriptutorial/blob/master/scripts/05_network_reconstruction_with_CARNIVAL.md


#Increase memory for CARNIVAL:
#library(unix)
#rlimit_as(9e+10)  #increases to ~90GB

#The whole script is for one selected brain region and one selected gender:

selected_brain_region="nAcc"
selected_gender="M"

#for exporting: 

condition="scz_vs_control"
sex_brain_region=paste(selected_gender, "_", selected_brain_region, sep="")
path=source_path
tag=paste0(sex_brain_region, "_", condition, "_dorothea_gene_symbol_file_input", sep="", collapse=NULL)


#Run CARNIVAL
#Import datasets: 
tf_activities <- read_csv(paste0(path, tag, "dorothea_TFActivity_CARNIVALinput.csv"))
PathwayActivity <- read_csv(paste0(path, tag, "PathwayActivity_CARNIVALinput.csv"))

#Before running CARNIVAL, we need to create or upload a scaffold network. 
#This will be “the map” that the ILP algorithm will follow to find the causal network. 
#We use Omnipath to obtain the signed and directed interactions from all the available resources. 

#CARNIVAL requires this information in a sif table (node1, interaction, node2) format, 
#therefore we use the consensus columns of direction (consensus_direction) 
#and sign (consensus_stimulation and consensus_inhibition) to extract it.

#The query returns 0/1 as logic status of being a stimulation or an inhibition reaction. 
#Thus, this output is reformulated as 1/-1 to indicate stimulation or inhibition, respectively. 
#We can keep either the interactions that are consistent, or both alternatives (e.g. A 1 B; A -1 B). 
#In this example, we keep the consistent ones.

omniR <- import_omnipath_interactions()

## Warning: 'import_Omnipath_Interactions' is deprecated.
## Use 'import_omnipath_interactions' instead.
## See help("Deprecated")

# signed and directed
omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
                                (consensus_stimulation == 1 | 
                                 consensus_inhibition == 1
                                 ))
  
# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which( omnipath_sd$consensus_stimulation == 0)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 1)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 0)] = 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
      dplyr::filter(consensus_stimulation==consensus_inhibition) %>%
      unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c('source', 'interaction', 'target')

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)


#save SIF
write_tsv(sif, paste0(path, tag,"omnipath_carnival_cplex_with_weights.tsv"))


#Transcription Factor and pathway activities for CARNIVAL

#We use the supplementary functions generateTFList.r and assignPROGENyScores.r 
#to shift the formats of tf_activities and PathwayActivity to the one required by CARNIVAL.

# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$TF
tf_activities_carnival$TF <- NULL
tfList = generateTFList(tf_activities_carnival, top=50, access_idx = 1)

# progeny for CARNIVAL
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
rownames(PathwayActivity_carnival) <- PathwayActivity_carnival$Pathway
PathwayActivity_carnival$Pathway <- NULL
progenylist = assignPROGENyScores(progeny = t(PathwayActivity_carnival), 
                                            progenyMembers = progenyMembers, 
                                            id = "gene", 
                                            access_idx = 1)

#Running CARNIVAL

#CARNIVAL has been developed to find the causal link between the activities of the transcription factors (TFs) and the ‘perturbed’ nodes. 
#In current version, v1.0.0, we have 3 main inputs that we have to provide:

    ### measObj: The TFs’ activities (like the ones we have obtained from DoRothEA)
    ### inputObj: The ‘perturbed’ nodes we want that CARNIVAL connects with the activity of TFs. 
    ### There are 3 ways of using it:

    ### 1 Give the name and sign of the selected nodes;
    ### 2 Give the name only, so the algorithm will select the sign that best fit the models,
    ### 3 Give NULL as value will create a “Perturbation” node that will try both signs for all ‘initial’ nodes of the given network ( netObj ).

    ### netObj: The network that will serve as map to connect the TFs’ activities ( measObj ) and the perturbed nodes ( inputObj )

#Although it is not required, a fourth object called weightObj can be also given. 
#This object gives values ranged from -1 to 1 for a set of nodes of the network. 
#The aim of weightObj is helping the solver to find optimal solutions faster.

#In the present example, we use assign as perturbation nodes all the “initial” nodes (option 2), 
#and as weightObj the PROGENy scores assigned to the most representative genes of the calculated pathways

# get initial nodes
iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

# run carnival
carnival_result = runCARNIVAL( inputObj= iniciators,
                               measObj = tfList$t, 
                               netObj = sif, 
                               weightObj = progenylist$score, 
                               solver = "cplex",
                               timelimit=7200,
                               mipGAP=0,
                               poolrelGAP=0 )


#resultsLpSolve <- runVanillaCarnival( perturbations = iniciators, 
#                                      measurements = tfList$t,
#                                      priorKnowledgeNetwork = sif, 
#                                      solver = "lpSolve",
#                                      timelimit=7200,
#                                      mipGAP=0,
#                                      poolrelGAP=0 )

#CARNIVAL gives a list of 4 elements:

    #weightedSIF: summary of all interactions found in all models
    #nodesAttributes: summary of all nodes and how many times are they found in the different models
    #sifAll: networks of all the models
    #attributesAll: node attributes of all models

#We can now visualise the network…

#transoform to data.frame
carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)

carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)

saveRDS(carnival_result, paste0(path, tag,"carnival_result_cplex_with_weights.rds"))

# visualization
visNet = carnival_visNet(evis = carnival_result$weightedSIF,
                         nvis = carnival_result$nodesAttributes)

## Graphical representation of sample
#visNet
visSave(visNet, file = paste0(path, tag,'carnival_visualization_visNetwork_cplex_with_weights.html'), selfcontained = TRUE)