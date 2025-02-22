library (mvMORPH)
library (phytools)
library (ape)
library (factoextra)
library (cluster)
library (webshot2)
library (rgl)
library (viridisLite)
library (corHMM)
library (dplyr)
library (purrr)
library (ggplot2)
library (tidyr)
library (ggpubr)
library (parallel)
library (reshape)
library (forcats)
library (plotly)
library (phyloglm)
library (ggh4x)
library (htmlwidgets)

source("./functions.R") #Various mildly-helpful functions

############################################################################
##### Reading RDS in case you don't feel like running all the analyses #####
############################################################################

folder_path <- "./RDS"
rds_files <- list.files(folder_path, pattern = "\\.RDS$", full.names = TRUE)
for(file_path in rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  data <- readRDS(file_path)
  assign(file_name, data, envir = .GlobalEnv)
}

#For reproducibility
set.seed(69420)
sampler <- sample(1:100, 1)
indices <- 1:100
treesampler <- 100

####################################
############ Read data #############
####################################

traits <- read.table("./Data/traits.txt", 
                     header = TRUE, 
                     sep = "\t")
wc <- read.table("./Data/worldclim.txt", 
                 header = TRUE, 
                 sep = "\t")
wc <- na.omit(wc)

####################################
############ Read trees ############
####################################

#Read and prepare MCMC trees
tree1 <- read.nexus("./Tree/MrBayes/Helio_trim.nex.run1.t")
tree2 <- read.nexus("./Tree/MrBayes/Helio_trim.nex.run2.t")
tree.merge <- c(tree1[round(length(tree1)/4):length(tree1)],
                tree2[round(length(tree2)/4):length(tree2)])
tree.rooted <- root.multiPhylo(tree.merge, 
                               "Chamira_circaeoides", 
                               resolve.root = TRUE)
tree.sample <- sample(tree.rooted, treesampler)
tree.sample <- lapply(tree.sample, ladderize)

contree <- read.nexus("./Tree/MrBayes/Helio_trim.nex.con.tre")
contree.rooted <- root(contree, "Chamira_circaeoides")
contree.rooted <- ladderize(contree.rooted)

saveRDS(tree.sample, "./RDS/tree.sample.RDS")

#################################################
############ Subset other dataframes ############
#################################################

anatomy <- traits[,c(1,8:16)]
anatomy.rnames <- anatomy
row.names(anatomy.rnames) <- anatomy.rnames[,1]
anatomy.rnames <- anatomy.rnames[,-1]
anatomy.scaled <- scale(anatomy.rnames)
anatomy.scaled <- anatomy.scaled[row.names(anatomy.scaled) %in% tree.sample[[1]]$tip.label,]

discrete <- traits[,1:3]

BIO <- wc[,c(1,38:ncol(wc))]
BIO <- aggregate(. ~ species, data = BIO, FUN = mean)
rownames(BIO) <- BIO$species
BIO <- BIO[,-1]
BIO.scaled <- scale(BIO)

#######################################
############ PCA - climate ############
#######################################

# This part of the script runs principal component analysis on bioclimatic
# variable, and prints some output statistics

commontips <- intersect(row.names(BIO.scaled), tree.sample[[1]]$tip.label)
BIO.scaled.PCA <- BIO.scaled[row.names(BIO.scaled) %in% commontips,]
tree.sample.PCA <- keep.tip.multiPhylo(tree.sample, commontips)

#Standard PCA

PCA <- prcomp(BIO.scaled, center = FALSE, scale. = FALSE) #already scaled
PCs <- PCA$x[,1:3]

#Variance by component

round(PCA$sdev^2 / sum(PCA$sdev^2)*100, 4)

#Make table with important components

threshold <- 0.3
PCA.vecs <- PCA$rotation[,1:3]
ifelse(abs(PCA.vecs) > threshold, PCA.vecs, NA)

#####################################################
################## ACE - discrete ###################
#####################################################

# This part of the script performs ancestral character estimation for lifespan
# and scondary woodiness using package corHMM. We are checking two possible
# transition rate matrices - ER (equal rates, assimuing single rate of transition 
#in both directions) and ARD (all rates different, assuming two transition rates)

############
#### ER ####
############

#Transition rate matrix

ER <- getStateMat4Dat(discrete[,1:2], model = "ER")$rate.mat

#Lifespan
corHMM.ER.ls <- list()
for(i in seq_along(tree.sample)) {
  corHMM.ER.ls[[i]] <- corHMM(phy = tree.sample[[i]], 
                              data = discrete[,c(1,3)], 
                              nstarts = 20, rate.cat = 1, 
                              rate.mat = ER, 
                              node.states = "marginal",
                              n.cores = 20)
}

#Woodiness
corHMM.ER.wd <- list()
for(i in seq_along(tree.sample)) {
  corHMM.ER.wd[[i]] <- corHMM(phy = tree.sample[[i]], 
                              data = discrete[,c(1,2)], 
                              nstarts = 20, rate.cat = 1, 
                              rate.mat = ER, 
                              node.states = "marginal",
                              n.cores = 20)
}


#############
#### ARD ####
#############

#Transition rate matrix

ARD <- getStateMat4Dat(discrete[,1:2], model = "ARD")$rate.mat

#Lifespan
corHMM.ARD.ls <- list()
for(i in seq_along(tree.sample)) {
  corHMM.ARD.ls[[i]] <- corHMM(phy = tree.sample[[i]], 
                               data = discrete[,c(1,3)], 
                               nstarts = 20, rate.cat = 1, 
                               rate.mat = ARD, 
                               node.states = "marginal",
                               n.cores = 20)
}

#Woodiness
corHMM.ARD.wd <- list()
for(i in seq_along(tree.sample)) {
  corHMM.ARD.wd[[i]] <- corHMM(phy = tree.sample[[i]], 
                               data = discrete[,c(1,2)], 
                               nstarts = 20, rate.cat = 1, 
                               rate.mat = ARD, 
                               node.states = "marginal",
                               n.cores = 20)
}

###############################################################
#### Compare AICc for best and generate plots with summary ####
###############################################################

#Lifespan
ls.summary <- data.frame(Tree = indices,
                         ER = AICc.vec(corHMM.ER.ls),
                         ARD = AICc.vec(corHMM.ARD.ls))

ls.summary.melt <- ls.summary %>%
  gather(key = "Model", value = "AICc", -Tree) 

ls.summary.melt <- ls.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, AICc)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ AICc,
           row_number() == 2 ~ AICc - lag(AICc),
           row_number() == 3 ~ AICc - lag(AICc, order_by = AICc, n = 1)
         )
  )

ls.summary.plot <- ggplot(ls.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  theme(legend.position="none") +
  scale_fill_manual(values = c("firebrick4", "springgreen4")) +
  coord_cartesian(ylim = c(min(ls.summary.melt$AICc), NA)) +
  # scale_y_reverse() +
  labs(title = "",
       x = "Tree", y = "AICc", fill = "Model") +
  theme_minimal() +
  scale_x_discrete(breaks = seq(from = 1, to = 100, by = 10)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

ggsave("corHMM.pdf", ls.summary.plot, device = "pdf", width = 2.5, height = 2, units = "in")

#Woodiness
wd.summary <- data.frame(Tree = indices,
                         ER = AICc.vec(corHMM.ER.wd),
                         ARD = AICc.vec(corHMM.ARD.wd))

wd.summary.melt <- wd.summary %>%
  gather(key = "Model", value = "AICc", -Tree) 

wd.summary.melt <- wd.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, AICc)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ AICc,
           row_number() == 2 ~ AICc - lag(AICc),
           row_number() == 3 ~ AICc - lag(AICc, order_by = AICc, n = 1)
         )
  )

wd.summary.plot <- ggplot(wd.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  theme(legend.position="none") +
  scale_fill_manual(values = c("firebrick4", "springgreen4")) +
  coord_cartesian(ylim = c(min(wd.summary.melt$AICc), NA)) +
  # scale_y_reverse() +
  labs(title = "",
       x = "Tree", y = "AICc", fill = "Model") +
  theme_minimal() +
  scale_x_discrete(breaks = seq(from = 1, to = 100, by = 10)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

ggsave("./Plots/corHMM_wood.pdf", wd.summary.plot, device = "pdf", width = 2.5, height = 2, units = "in")

#######################
#### FOR CONSENSUS ####
#######################

#This part conducts reocnstruction for best model on the consensus tree and plots the output

#Lifespan
corHMM.ER.ls.plot <- corHMM(phy = contree.rooted, 
                            data = discrete[,c(1,3)], 
                            nstarts = 20, rate.cat = 1, 
                            rate.mat = ER, 
                            node.states = "marginal",
                            n.cores = 20)

sorted.ls <- discrete[match(corHMM.ER.ls.plot$phy$tip.label, discrete$species),] 
tip.colors.ls <- ifelse(sorted.ls$lifespan == "?", "grey30", 
                        ifelse(sorted.ls$lifespan == 0, "goldenrod1", "darkorchid3"))

pdf("./Plots/marginal_lifespan", height = 8.66, width = 5.00)
par(lwd = 0.01)
plotRECONJB(corHMM.ER.ls.plot$phy, 
            corHMM.ER.ls.plot$states, 
            show.tip.label = TRUE,
            tipcolors = tip.colors.ls,
            pie.cex = 0.5, 
            piecolors = c("goldenrod1", "darkorchid3"), 
            height = 20, width = 11, 
            edge.width = 0.5)
axisPhylo(cex = 0.3)
par(mar=c(5,3,2,2)+0.1)
dev.off()

#Woodiness
corHMM.ARD.wd.plot <- corHMM(phy = contree.rooted, 
                             data = discrete[,c(1,2)], 
                             nstarts = 20, rate.cat = 1, 
                             rate.mat = ARD, 
                             node.states = "marginal",
                             n.cores = 20)

sorted.wd <- discrete[match(corHMM.ARD.wd.plot$phy$tip.label, discrete$species),] 
tip.colors.wd <- ifelse(sorted.wd$wood == "?", "grey30", 
                        ifelse(sorted.wd$wood == 0, "salmon", "dodgerblue"))

pdf("./Plots/marginal_wood.pdf", height = 8.66, width = 3.50)
par(lwd = 0.01)
plotRECONJB(corHMM.ARD.wd.plot$phy, 
            corHMM.ARD.wd.plot$states, 
            show.tip.label = TRUE,
            tipcolors = tip.colors.wd,
            pie.cex = 0.5,
            cex = 0.4,
            piecolors = c("salmon", "dodgerblue"), 
            height = 20, width = 11, 
            edge.width = 0.5)
axisPhylo(cex = 0.3)
par(mar=c(5,3,2,2)+0.1)
dev.off()

#############################################################
#############################################################
################## Hypotheses testing <3 ####################
#############################################################
#############################################################

#Now we can finally move to the main part of the analyses!

##########################
##### Create SIMMAPs #####
##########################

#We really need stochastic maps only for the lifespan, but my junky
#code relies on some outputs for wood anatomy as well. Sorry!

simmaps.ls <- list() 
for(i in seq_along(tree.sample)) {
  simmaps.ls[[i]] <- makeSimmap(tree.sample[[i]],
                                data = discrete[,c(1,3)],
                                collapse = FALSE,
                                nSim = 1,
                                model = corHMM.ER.ls[[i]]$solution,
                                rate.cat = 1,
                                nCores = 20)
}

saveRDS(simmaps.ls, "./RDS/simmaps.ls.RDS")

simmaps.wd <- list() 
for(i in seq_along(tree.sample)) {
  simmaps.wd[[i]] <- makeSimmap(tree.sample[[i]],
                                data = discrete[,c(1,2)],
                                collapse = FALSE,
                                nSim = 1,
                                model = corHMM.ARD.wd[[i]]$solution,
                                rate.cat = 1,
                                nCores = 20)
}

saveRDS(simmaps.wd, "./RDS/simmaps.wd.RDS")

######################################
#### Woodiness/Lifespan ~ Climate ####
######################################

# In this part of the script, we are using phylogenetic logistic regression to 
# test for the association between climate (PCs reconstructed earlier) and lifespan
# or wood anatomy

#Extract lifespan adn wood data from the original dataset
lifespan <- data.frame(
  row.names = traits$species,
  lifespan = traits [,3])
lifespan <- subset (lifespan, lifespan != "?")

wood <- data.frame(
  row.names = traits$species,
  wood = traits [,2])
wood <- subset (wood, wood != "?")

#Prepare dataset with PCs
PCs.anatols <- PCs [row.names(PCs) %in% row.names(lifespan),]
rowdiff <- setdiff(simmaps.ls[[1]][[1]]$tip.label, row.names(PCs.anatols))
new.rows <- matrix (NA, nrow = length(rowdiff), ncol = ncol(PCs.anatols)) 
row.names(new.rows) <- rowdiff
PCs.anatols.fin <- na.omit(as.data.frame(rbind (PCs.anatols, new.rows)))
PCs.anatols.fin$lifespan <- lifespan[rownames(PCs.anatols.fin), "lifespan"]

PCs.anatowd <- PCs [row.names(PCs) %in% row.names(wood),]
rowdiff <- setdiff(simmaps.wd[[1]][[1]]$tip.label, row.names(PCs.anatowd))
new.rows <- matrix (NA, nrow = length(rowdiff), ncol = ncol(PCs.anatowd)) 
row.names(new.rows) <- rowdiff
PCs.anatowd.fin <- na.omit(as.data.frame(rbind (PCs.anatowd, new.rows)))
PCs.anatowd.fin$wood <- wood[rownames(PCs.anatowd.fin), "wood"]
PCs.anatowd.fin <-  PCs.anatowd.fin[row.names(PCs.anatowd.fin) %in% traits[traits$wood == 1,]$species,]

#Subset tip labels to match the datasets
tree.sample.lifespanoclim <- keep.tip.multiPhylo(tree.sample.PCA, row.names(PCs.anatols.fin))
tree.sample.woodoclim <- keep.tip.multiPhylo(tree.sample.PCA, row.names(PCs.anatowd.fin))

#Merge trait data with PCs
PCs.anatols.fin$lifespan <- lifespan[rownames(PCs.anatols.fin), "lifespan"]
PCs.anatowd.fin$wood <- wood[rownames(PCs.anatowd.fin), "wood"]

#Run phyloglm
lifespanoclim <- mclapply(indices, function(i) {
  phyloglm(lifespan ~ PC1 + PC2 + PC3,
           data = PCs.anatols.fin,
           tree.sample.lifespanoclim[[i]],
           method = "logistic_MPLE",
           boot = 1000
  )
}, mc.cores = 20)

woodoclim <- mclapply(indices, function(i) {
  phyloglm(wood ~ PC1 + PC2 + PC3,
           data = PCs.anatowd.fin,
           tree.sample.woodoclim[[i]],
           method = "logistic_MPLE",
           boot = 1000
  )
}, mc.cores = 20)

#Plot p-values for regression

#Lifespan
lifespanoclim.pval <- data.frame(Tree = indices, 
                                 p_value_PC1 = pvalvec_phylolm(lifespanoclim, 2),
                                 p_value_PC2 = pvalvec_phylolm(lifespanoclim, 3),
                                 p_value_PC3 = pvalvec_phylolm(lifespanoclim, 4))
lifespanoclim.pval_long <- lifespanoclim.pval %>%
  pivot_longer(cols = starts_with("p_value_PC"), 
               names_to = "Predictor", 
               values_to = "p_value") %>%
  mutate(Predictor = recode(Predictor, 
                            "p_value_PC1" = "PC1",
                            "p_value_PC2" = "PC2",
                            "p_value_PC3" = "PC3"))
lifespanoclim.pval_long <- lifespanoclim.pval_long %>%
  group_by(Tree) %>%
  mutate(x_position = as.numeric(Tree) + (as.numeric(factor(Predictor)) - 2) * 0.2)  # Offset for predictors
max_p_value <- max(lifespanoclim.pval_long$p_value, na.rm = TRUE)

#Woodiness
woodoclim.pval <- data.frame(Tree = indices, 
                             p_value_PC1 = pvalvec_phylolm(woodoclim, 2),
                             p_value_PC2 = pvalvec_phylolm(woodoclim, 3),
                             p_value_PC3 = pvalvec_phylolm(woodoclim, 4))
woodoclim.pval_long <- woodoclim.pval %>%
  pivot_longer(cols = starts_with("p_value_PC"), 
               names_to = "Predictor", 
               values_to = "p_value") %>%
  mutate(Predictor = recode(Predictor, 
                            "p_value_PC1" = "PC1",
                            "p_value_PC2" = "PC2",
                            "p_value_PC3" = "PC3"))
woodoclim.pval_long <- woodoclim.pval_long %>%
  group_by(Tree) %>%
  mutate(x_position = as.numeric(Tree) + (as.numeric(factor(Predictor)) - 2) * 0.2)  # Offset for predictors
max_p_value <- max(woodoclim.pval_long$p_value, na.rm = TRUE)

#Generate plots and save as pdfs
lifespanoclim.pval.plot <- ggplot(lifespanoclim.pval_long, aes(x = x_position, y = p_value, fill = Predictor)) +
  geom_bar(stat = "identity", width = 0.2, position = position_identity()) +  # Keep bars within a tree close
  geom_hline(yintercept = 0.05, color = "black", linetype = "solid") +  # Significance line
  labs(title = "", x = "Tree", y = "p-value") +
  coord_cartesian(ylim = c(0, max_p_value)) +  # Dynamically adjust y axis range
  scale_fill_manual(values = c("PC1" = "firebrick4", "PC2" = "darkorange", "PC3" = "springgreen4")) +  # Custom colors
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_continuous(breaks = seq(1, 100, by = 10),  # Show every 10th tree
                     labels = seq(1, 100, by = 10),   # Tree labels
                     expand = c(0, 0)) +  # Ensure axes are tight against the bars
  scale_y_continuous(expand = c(0, 0))  # Ensure y-axis is tight against the bars

ggsave("./Plots/lifespanoclim_pval.pdf", 
       lifespanoclim.pval.plot, 
       device = "pdf", 
       width = 6, 
       height = 1.5, 
       units = "in")

woodoclim.pval.plot <- ggplot(woodoclim.pval_long, aes(x = x_position, y = p_value, fill = Predictor)) +
  geom_bar(stat = "identity", width = 0.2, position = position_identity()) +  # Keep bars within a tree close
  geom_hline(yintercept = 0.05, color = "black", linetype = "solid") +  # Significance line
  labs(title = "", x = "Tree", y = "p-value") +
  coord_cartesian(ylim = c(0, max_p_value)) +  # Dynamically adjust y axis range
  scale_fill_manual(values = c("PC1" = "firebrick4", "PC2" = "darkorange", "PC3" = "springgreen4")) +  # Custom colors
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_continuous(breaks = seq(1, 100, by = 10),  # Show every 10th tree
                     labels = seq(1, 100, by = 10),   # Tree labels
                     expand = c(0, 0)) +  # Ensure axes are tight against the bars
  scale_y_continuous(expand = c(0, 0))  # Ensure y-axis is tight against the bars

ggsave("./Plots/woodoclim_pval.pdf", 
       woodoclim.pval.plot, 
       device = "pdf", 
       width = 6, 
       height = 1.5, 
       units = "in")

#Plot the output regression curves

#Lifespan
plot_logistic_curve_ggplot_ls <- function(predictor_name, predictor_values, response_values, models_list, predictor) {
  
  predictor_values <- as.numeric(predictor_values)
  response_values <- as.numeric(response_values)
  
  plot_data <- data.frame(Predictor = predictor_values, Response = response_values)
  
  curve_data <- data.frame()
  
  extended_range <- range(predictor_values)
  x_values <- seq(extended_range[1] - 1, extended_range[2] + 1, length.out = 100)  # Extended range by 1
  
  for (i in 1:length(models_list)) {
    fit <- models_list[[i]]
    cc <- coef(fit)
    
        if (predictor == "PC1") {
      coef_index <- 2
    } else if (predictor == "PC2") {
      coef_index <- 3 
    } else if (predictor == "PC3") {
      coef_index <- 4
    } else {
      next  # Skip if an unexpected predictor is passed
    }
    
    # Apply the logistic function with the correct coefficient for the predictor
    y_values <- plogis(cc[1] + cc[coef_index] * x_values)
    curve_data <- rbind(curve_data, data.frame(X = x_values, Y = y_values, Tree = i))
  }
  
  curve_color <- switch(predictor,
                        "PC1" = "firebrick4",
                        "PC2" = "darkorange",
                        "PC3" = "springgreen4")
  
  p <- ggplot(plot_data, aes(x = Predictor, y = Response, color = as.factor(Response))) +
    geom_jitter(height = 0.02, width = 0, alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("0" = "goldenrod1", "1" = "darkorchid3")) +  # Set colors for annuals and perennials
    scale_y_continuous(breaks = c(0, 1), limits = c(-0.02, 1.02)) +
    labs(x = predictor_name, y = "Lifespan") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Add all logistic regression curves using geom_line with the specified color
  p <- p + geom_line(data = curve_data, aes(x = X, y = Y, group = Tree), color = curve_color, alpha = 0.2)
  return(p)
}


#Woodiness
plot_logistic_curve_ggplot_wd <- function(predictor_name, predictor_values, response_values, models_list, predictor) {
  
  predictor_values <- as.numeric(predictor_values)
  response_values <- as.numeric(response_values)
  
  plot_data <- data.frame(Predictor = predictor_values, Response = response_values)
  
  curve_data <- data.frame()
  
  # Extend the range of x-values beyond the observed range
  extended_range <- range(predictor_values)
  x_values <- seq(extended_range[1] - 1, extended_range[2] + 1, length.out = 100)  # Extended range by 1
  
  for (i in 1:length(models_list)) {
    fit <- models_list[[i]]
    cc <- coef(fit)
    
    # Check for the correct coefficient position based on the predictor
    if (predictor == "PC1") {
      coef_index <- 2  # PC1
    } else if (predictor == "PC2") {
      coef_index <- 3  # PC2
    } else if (predictor == "PC3") {
      coef_index <- 4  # PC3
    } else {
      next  # Skip if an unexpected predictor is passed
    }
    
    # Apply the logistic function with the correct coefficient for the predictor
    y_values <- plogis(cc[1] + cc[coef_index] * x_values)
    curve_data <- rbind(curve_data, data.frame(X = x_values, Y = y_values, Tree = i))
  }
  
  curve_color <- switch(predictor,
                        "PC1" = "firebrick4",
                        "PC2" = "darkorange",
                        "PC3" = "springgreen4")
  
  p <- ggplot(plot_data, aes(x = Predictor, y = Response, color = as.factor(Response))) +
    geom_jitter(height = 0.02, width = 0, alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("0" = "salmon", "1" = "dodgerblue")) +  # Set colors for annuals and perennials
    scale_y_continuous(breaks = c(0, 1), limits = c(-0.02, 1.02)) +  # Show only 0 and 1 on the y-axis
    labs(x = predictor_name, y = "Secondary woodiness") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Add all logistic regression curves using geom_line with the specified color
  p <- p + geom_line(data = curve_data, aes(x = X, y = Y, group = Tree), color = curve_color, alpha = 0.2)
  
  return(p)
}

#Generate plots and save outputs as PDFs
PC1_values_ls <- PCs.anatols.fin$PC1
PC2_values_ls <- PCs.anatols.fin$PC2
PC3_values_ls <- PCs.anatols.fin$PC3
lifespan_values <- PCs.anatols.fin$lifespan 

p1_ls <- plot_logistic_curve_ggplot_ls("PC1", PC1_values_ls, lifespan_values, lifespanoclim, "PC1")
p2_ls <- plot_logistic_curve_ggplot_ls("PC2", PC2_values_ls, lifespan_values, lifespanoclim, "PC2")
p3_ls <- plot_logistic_curve_ggplot_ls("PC3", PC3_values_ls, lifespan_values, lifespanoclim, "PC3")

PC1_values_wd <- PCs.anatowd.fin$PC1
PC2_values_wd <- PCs.anatowd.fin$PC2
PC3_values_wd <- PCs.anatowd.fin$PC3
wood_values <- PCs.anatowd.fin$wood  # Binary response (0 = annual, 1 = perennial)

p1_wd <- plot_logistic_curve_ggplot_wd("PC1", PC1_values_wd, wood_values, woodoclim, "PC1")
p2_wd <- plot_logistic_curve_ggplot_wd("PC2", PC2_values_wd, wood_values, woodoclim, "PC2")
p3_wd <- plot_logistic_curve_ggplot_wd("PC3", PC3_values_wd, wood_values, woodoclim, "PC3")

ggsave("./Plots/regression_PC1.pdf", 
       p1_ls, 
       device = "pdf", 
       width = 2, 
       height = 2.5, 
       units = "in")

ggsave("./Plots/regression_PC2.pdf", 
       p2_ls, 
       device = "pdf", 
       width = 2, 
       height = 2.5, 
       units = "in")

ggsave("./Plots/regression_PC3.pdf", 
       p3_ls, 
       device = "pdf", 
       width = 2, 
       height = 2.5, 
       units = "in")

ggsave("./Plots/regression_PC1_wood.pdf", 
       p1_wd, 
       device = "pdf", 
       width = 2, 
       height = 2.5, 
       units = "in")

ggsave("./Plots/regression_PC2_wood.pdf", 
       p2_wd, 
       device = "pdf", 
       width = 2, 
       height = 2.5, 
       units = "in")

ggsave("./Plots/regression_PC3_wood.pdf", 
       p3_wd, 
       device = "pdf", 
       width = 2, 
       height = 2.5, 
       units = "in")

############################
#### Anatomy ~ lifespan ####
############################

# This part of the script fits four models BM1, BMM, OU1 and OUM to continuous wood anatomical traits
# using previously reconstructed lifespan simmap to specify selective regimes. It takes a while to finish!

anatols.BM1 <- mclapply(indices, function(i) {
  mvBM(simmaps.ls[[i]][[1]], 
       anatomy.scaled, 
       model = "BM1", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(anatols.BM1, "./RDS/anatols.BM1.RDS")

anatols.BMM <- mclapply(indices, function(i) {
  mvBM(simmaps.ls[[i]][[1]], 
       anatomy.scaled, 
       model = "BMM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(anatols.BMM, "./RDS/anatols.BMM.RDS")

anatols.OU1 <- mclapply(indices, function(i) {
  mvOU(simmaps.ls[[i]][[1]], 
       anatomy.scaled, 
       model = "OU1", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(anatols.OU1, "./RDS/anatols.OU1.RDS")

anatols.OUM <- mclapply(indices, function(i) {
  mvOU(simmaps.ls[[i]][[1]], 
       anatomy.scaled, 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(anatols.OUM, "./RDS/anatols.OUM.RDS")

AICc.compare(anatols.BM1)
AICc.compare(anatols.BMM)
AICc.compare(anatols.OU1)
AICc.compare(anatols.OUM)

#Now, we can create the barplots showing model fit across 100 trees
anatols.summary <- data.frame(Tree = indices,
                              BM1 = AICc.vec(anatols.BM1),
                              BMM = AICc.vec(anatols.BMM),
                              OU1 = AICc.vec(anatols.OU1),
                              OUM = AICc.vec(anatols.OUM))

anatols.summary.melt <- anatols.summary %>%
  gather(key = "Model", value = "AICc", -Tree) 

anatols.summary.melt <- anatols.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, AICc)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ AICc,
           row_number() == 2 ~ AICc - lag(AICc),
           row_number() == 3 ~ AICc - lag(AICc, order_by = AICc, n = 1),
           row_number() == 4 ~ AICc - lag(AICc, order_by = AICc, n = 2)
         )
  )

anatols.summary.plot <- ggplot(anatols.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_discrete(breaks = levels(factor(anatols.summary.melt$Tree))[seq(1, length(levels(factor(anatols.summary.melt$Tree))), 10)]) +
  scale_fill_manual(values = c("firebrick4", "indianred", "springgreen4", "palegreen")) +
  coord_cartesian(ylim = c(min(anatols.summary.melt$AICc), NA)) +
  theme(legend.position="none") +
  # scale_y_reverse() +
  labs(title = element_blank(),
       x = "Tree", y = "AICc", fill = "Model")

ggsave("./Plots/anatols.pdf", 
       anatols.summary.plot, 
       device = "pdf", 
       width = 1.5, 
       height = 6, 
       units = "in")

# Finally, we can create density plots for evolutionary optima. Since model fitting
# was kind of inonclusive wherer OU1 or OUM fits better, we generate plots for both scenarios

#Revert scaling

anatomean <- colMeans(na.omit(anatomy.rnames))
anatoscale <- sqrt(diag(var(na.omit(anatomy.rnames))))

#OUM
anatols.OUM.theta <- theta.compare(anatols.OUM, anatoscale, anatomean)

anatols.OUM.theta.bind <- theta_converter(anatols.OUM.theta, "theta", "lifespan")

anatols.OUM.theta.bind <- anatols.OUM.theta.bind %>%
  mutate(column = factor(lifespan), 
         panel_id = factor(panel_id),
         list_id = factor(list_id))

anatols.OUM.theta.bind$lifespan <- gsub("0", "annuals", anatols.OUM.theta.bind$lifespan)
anatols.OUM.theta.bind$lifespan <- gsub("1", "perennials", anatols.OUM.theta.bind$lifespan)

anatols.OUM.theta.plot <- ggplot(anatols.OUM.theta.bind, aes(x = value, fill = interaction(list_id, lifespan))) + 
  geom_density(alpha = 0.6) + # Add some transparency
  facet_wrap(~row_name, scales = "free", ncol = 1) + # Create a panel for each row, adjust ncol as needed
  labs(title = element_blank(), 
       x = "", 
       y = "Density") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_manual(values = c("goldenrod1", "darkorchid1")) +
  theme(legend.position = "none") +
  labs(fill = "")

ggsave("./Plots/anatols_OUM_theta.pdf", 
       anatols.theta.plot, 
       device = "pdf", 
       width = 3, 
       height = 6, 
       units = "in")

#OU1
anatols.OU1.theta <- theta.compare(anatols.OU1, anatoscale, anatomean)

anatols.OU1.theta.bind <- theta_converter_0(anatols.OU1.theta, "theta", "OU1")

anatols.OU1.theta.bind <- anatols.OU1.theta.bind %>%
  mutate(column = factor(OU1), 
         panel_id = factor(panel_id),
         list_id = factor(list_id))

range_oum <- anatols.OUM.theta.bind %>%
  group_by(row_name) %>%
  summarise(min_value_oum = min(value), max_value_oum = max(value), .groups = "drop")

# Join the OUM-based limits to the OU1 dataset
anatols.OU1.theta.bind <- anatols.OU1.theta.bind %>%
  left_join(range_oum, by = "row_name")

scales_x <- lapply(1:nrow(range_oum), function(i) {
  scale_x_continuous(limits = c(range_oum$min_value_oum[i], range_oum$max_value_oum[i]))
})

# Create the plot for OU1 data with facet-specific x-axis limits
anatols.OU1.theta.plot <- ggplot(anatols.OU1.theta.bind, aes(x = value, fill = interaction(list_id, OU1))) + 
  geom_density(alpha = 0.6) + 
  facet_wrap(~row_name, scales = "free", ncol = 1) +  # Enable free x and y-scales
  labs(title = element_blank(), x = "", y = "Density") +
  theme_minimal() +
  theme(axis.line = element_line(color = 'black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_manual(values = c("chocolate3")) +
  theme(legend.position = "none") +
  facetted_pos_scales(x = scales_x, y = NULL)  # Apply custom x-axis scales, leave y free

# Save the plot
ggsave("./Plots/anatols_OU1_theta_facetted_scales.pdf", 
       anatols.OU1.theta.plot, 
       device = "pdf", 
       width = 3, 
       height = 6, 
       units = "in")

######################################
#### Anatomy ~ Lifespan + Climate ####
######################################
 
# In the last step we can fit phylogenetic regression models to test association between
# anatomy and climate using lifespan as a covariate. We are testing both OU and BM models.

#Fit GLS models - flowering season means

anatoclim <- merge(anatomy.scaled, PCs, by = "row.names")
colnames(anatoclim)[1] <- "ID"  # Rename the first column to avoid duplication

anatoclim <- merge(anatoclim, lifespan, by.x = "ID", by.y = "row.names")
row.names(anatoclim) <- anatoclim$ID
anatoclim <- anatoclim[, -1]
anatoclim <- na.omit(anatoclim)

#anatoclim [row.names (anatoclim) %in% traits[traits$lifespan == "0",]$species,]

anatoclim.glm <- list(anatomy = as.matrix(anatoclim[,1:9]), 
                      PC1 = as.matrix(anatoclim[,10]),
                      PC2 = as.matrix(anatoclim[,11]),
                      PC3 = as.matrix(anatoclim[,12]),
                      lifespan = as.matrix(anatoclim[,13]))                     

tree.sample.reduced.anatoclim <- keep.tip.multiPhylo(tree.sample, row.names(anatoclim))

anatoclim.fit.BM <- mclapply(indices, function(i) {
  mvgls(anatomy~PC1+PC2+PC3+PC1:lifespan+PC2:lifespan+PC3:lifespan+lifespan,
        data = anatoclim.glm,
        model = "BM",
        error = TRUE,
        tree.sample.reduced.anatoclim[[i]], 
        penalty = "RidgeArch")
}, mc.cores = 16)

saveRDS(anatoclim.fit.BM, "./RDS/anatoclim.fit.BM.RDS")

anatoclim.fit.OU <- mclapply(indices, function(i) {
  mvgls(anatomy~PC1+PC2+PC3+PC1:lifespan+PC2:lifespan+PC3:lifespan+lifespan, 
        data = anatoclim.glm,
        model = "OU",
        error = TRUE,
        tree.sample.reduced.anatoclim[[i]], 
        penalty = "RidgeArch")
}, mc.cores = 16)

saveRDS(anatoclim.fit.OU, "./RDS/anatoclim.fit.OU.RDS")

#Now we test model fit for both OU and BM using EIC

anatoclim.fit.BM.EIC <- list()
for(i in seq_along(anatoclim.fit.BM)){
  anatoclim.fit.BM.EIC[[i]] <- EIC(anatoclim.fit.BM[[i]])
}

anatoclim.fit.OU.EIC <- list()
for(i in seq_along(anatoclim.fit.OU)){
  anatoclim.fit.OU.EIC[[i]] <- EIC(anatoclim.fit.OU[[i]])
}

anatoclim.EIC.summary <- data.frame(Tree = indices,
                                    BM = EICvec(anatoclim.fit.BM.EIC),
                                    OU = EICvec(anatoclim.fit.OU.EIC))

anatoclim.EIC.summary.melt <- anatoclim.EIC.summary %>%
  gather(key = "Model", value = "anatoclim.EIC", -Tree) 

anatoclim.EIC.summary.melt  <- anatoclim.EIC.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, anatoclim.EIC)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ anatoclim.EIC,
           row_number() == 2 ~ anatoclim.EIC - lag(anatoclim.EIC),
           row_number() == 3 ~ anatoclim.EIC - lag(anatoclim.EIC, order_by = anatoclim.EIC, n = 1)
         )
  )

anatoclim.EIC.summary.plot <- ggplot(anatoclim.EIC.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_discrete(breaks = levels(factor(anatols.summary.melt$Tree))[seq(1, length(levels(factor(anatols.summary.melt$Tree))), 10)]) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("firebrick4", "springgreen4")) +
  coord_cartesian(ylim = c(min(anatoclim.EIC.summary.melt$anatoclim.EIC), NA)) +
  # scale_y_reverse() +
  labs(title = element_blank(),
       x = "Tree", y = "EIC", fill = "Model")

ggsave("./Plots/anatoclim.EIC.summary.plot.pdf", 
       anatoclim.EIC.summary.plot, 
       device = "pdf", 
       width = 3, 
       height = 6, 
       units = "in")

#MANOVA to test the signifficance

anatoclim.manova.BM <- mclapply(indices, function(i) {
  manova.gls(anatoclim.fit.BM[[i]], test = "Pillai", type = "III", nperm = 1000)
}, mc.cores = 16)

saveRDS(anatoclim.manova.BM, "./RDS/anatoclim.manova.BM.RDS")

anatoclim.manova.OU <- mclapply(indices, function(i) {
  manova.gls(anatoclim.fit.OU[[i]], test = "Pillai", type = "III", nperm = 1000)
}, mc.cores = 16)

saveRDS(anatoclim.manova.OU, "./RDS/anatoclim.manova.OU.RDS")

##############################################################
##############################################################
################## Various summary tables ####################
##############################################################
##############################################################

################################
#### Mean and SD for optima ####
################################

anatols.OUM.theta.meansd <- anatols.OUM.theta.bind %>%
  group_by(row_name, lifespan) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

saveRDS(anatols.OUM.theta.meansd, "./RDS/anatols.OUM.theta.meansd.RDS")

anatols.OU1.theta.meansd <- anatols.OU1.theta.bind %>%
  group_by(row_name, OU1) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

saveRDS(anatols.theta.meansd, "./RDS/anatols.OU1.theta.meansd.RDS")

###################
#### Halflives ####
###################

anatols.OUM.hl <- cbind.data.frame (meand = rowMeans(sapply(anatols.OUM, halflife)),
                                    sd = apply(sapply(anatols.OUM, halflife), 1, sd))
rownames(anatols.OUM.hl) <- colnames(anatols.OUM[[4]]$theta)

anatols.OU1.hl <- cbind.data.frame (meand = rowMeans(sapply(anatols.OU1, halflife)),
                                    sd = apply(sapply(anatols.OU1, halflife), 1, sd))
rownames(anatols.OU1.hl) <- colnames(anatols.OU1[[4]]$theta)

halflife.fin <-list (anatols.OU1.hl, anatols.OUM.hl)

saveRDS(halflife.fin, "./RDS/halflife.fin.RDS")

###############################################
#### Anatomy ~ Climate + Lifespan, summary ####
###############################################

############
#### BM ####
############

coefficients_list_BM <- lapply(indices, function(i) {
  anatoclim.fit.BM[[i]]$coefficients[2:8,]  # Extract rows for PC1, PC2, PC3 and their interactions with lifespan
})

pvalues_list_BM <- lapply(1:100, function(i) {
  anatoclim.manova.BM[[i]]$pvalue[2:8] 
})

# Step 2: Combine all coefficients into a 3D array (7 predictors x 9 traits x 100 trees)
coefficients_array_BM <- array(unlist(coefficients_list_BM), dim = c(7, 9, 100))
coefficients_mean_BM <- apply(coefficients_array_BM, c(1, 2), mean, na.rm = TRUE)
coefficients_sd_BM <- apply(coefficients_array_BM, c(1, 2), sd, na.rm = TRUE)
pvalues_matrix_BM <- do.call(rbind, pvalues_list_BM)

# Step 5: Average p-values and calculate SDs across 100 trees
pvalues_mean_BM <- apply(pvalues_matrix_BM, 2, mean, na.rm = TRUE)
pvalues_sd_BM <- apply(pvalues_matrix_BM, 2, sd, na.rm = TRUE)

# Step 6: Prepare the final summary table with mean ± SD for p-values and coefficients
predictors <- c("PC1:annual", "PC2:annual", "PC3:annual", "lifespan", "PC1:perennial", "PC2:perennial", "PC3:perennial")

# Create a data frame for the final table
summary_table_BM <- data.frame(
  Predictor = predictors,
  P_Value = paste0(round(pvalues_mean_BM, 3), " ± ", round(pvalues_sd_BM, 3)),
  t(sapply(1:7, function(i) {
    sapply(1:9, function(j) {
      paste0(round(coefficients_mean_BM[i, j], 3), " ± ", round(coefficients_sd_BM[i, j], 3))
    })
  }))
)

# Set column names for the traits
colnames(summary_table_BM)[-1] <- c("p-value", colnames(anatoclim[,1:9]))

saveRDS(summary_table_BM, "./RDS/summary_table_BM.RDS")

write.table(t(summary_table_BM), 
            "./Plots/summary_table_BM.tsv", 
            sep = "\t")

############
#### OU ####
############

coefficients_list_OU <- lapply(indices, function(i) {
  anatoclim.fit.OU[[i]]$coefficients[2:8,]  # Extract rows for PC1, PC2, PC3 and their interactions with lifespan
})

pvalues_list_OU <- lapply(1:100, function(i) {
  anatoclim.manova.OU[[i]]$pvalue[2:8] 
})

# Step 2: Combine all coefficients into a 3D array (7 predictors x 9 traits x 100 trees)
coefficients_array_OU <- array(unlist(coefficients_list_OU), dim = c(7, 9, 100))
coefficients_mean_OU <- apply(coefficients_array_OU, c(1, 2), mean, na.rm = TRUE)
coefficients_sd_OU <- apply(coefficients_array_OU, c(1, 2), sd, na.rm = TRUE)
pvalues_matrix_OU <- do.call(rbind, pvalues_list_OU)

# Step 5: Average p-values and calculate SDs across 100 trees
pvalues_mean_OU <- apply(pvalues_matrix_OU, 2, mean, na.rm = TRUE)
pvalues_sd_OU <- apply(pvalues_matrix_OU, 2, sd, na.rm = TRUE)

# Step 6: Prepare the final summary table with mean ± SD for p-values and coefficients
predictors <- c("PC1:annual", "PC2:annual", "PC3:annual", "lifespan", "PC1:perennial", "PC2:perennial", "PC3:perennial")

# Create a data frame for the final table
summary_table_OU <- data.frame(
  Predictor = predictors,
  P_Value = paste0(round(pvalues_mean_OU, 3), " ± ", round(pvalues_sd_OU, 3)),
  t(sapply(1:7, function(i) {
    sapply(1:9, function(j) {
      paste0(round(coefficients_mean_OU[i, j], 3), " ± ", round(coefficients_sd_OU[i, j], 3))
    })
  }))
)

# Set column names for the traits
colnames(summary_table_OU)[-1] <- c("p-value", colnames(anatoclim[,1:9]))

saveRDS(summary_table_OU, "./RDS/summary_table_OU.RDS")

write.table(t(summary_table_OU), 
            "./Plots/summary_table_OU.tsv", 
            sep = "\t")

#####################################################
#####################################################
################## Supplementary ####################
#####################################################
#####################################################

##################
#### PCA plot ####
##################

#Define clades

clade_A <- extract.clade(contree.rooted, 96)$tip.label
clade_B <- extract.clade(contree.rooted, 156)$tip.label
clade_C <- extract.clade(contree.rooted, 131)$tip.label
clade_D <- extract.clade(contree.rooted, 91)$tip.label
clade_E <- "Heliophila_alpina"

clades <- data.frame(
  row.names = c(clade_A, clade_B, clade_C, clade_D, clade_E),
  clade = c (rep("A", length(clade_A)),
             rep("B", length(clade_B)),
             rep("C", length(clade_C)),
             rep("D", length(clade_D)),
             rep("E", length(clade_E)))
)

PCs.plot <- merge(PCs, clades, by = "row.names")
row.names(PCs.plot) <- PCs.plot$Row.names
PCs.plot <- PCs.plot [,-1]
PCs.plot <- merge(PCs.plot, lifespan, by = "row.names")
row.names(PCs.plot) <- PCs.plot$Row.names
PCs.plot <- PCs.plot [,-1]

#PLot

# Create the 3D scatter plot

PCA_bioclim_plot <- plot_ly(PCs.plot, 
                            x = ~PC1, 
                            y = ~PC2, 
                            z = ~PC3, 
                            type = 'scatter3d', 
                            mode = 'markers',
                            color = ~factor(lifespan), 
                            colors = c("goldenrod1", "darkorchid1"), 
                            symbol = ~factor(clade),
                            symbols = c('circle', 'square', 'diamond', 'cross', 'x'),
                            marker = list(size = 4))


# Customize the layout

PCA_bioclim_plot <- PCA_bioclim_plot %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                                             yaxis = list(title = 'PC2'),
                                                             zaxis = list(title = 'PC3')),
                                                legend = list(title = list(text = "Lifespan & Clade")))

saveWidget(PCA_bioclim_plot, "./Plots/PCA_bioclim_plot.html")


######################
#### PCA treeplot ####
######################

contree.rooted.PCA <- keep.tip(contree.rooted, row.names(PCs.plot))

# Ensure species in PCs.plot match tip labels in the tree
PCs.plot$species <- rownames(PCs.plot)
species_order <- contree.rooted.PCA$tip.label
PCs.plot <- PCs.plot[PCs.plot$species %in% species_order, ]
PCs.plot <- PCs.plot[match(species_order, PCs.plot$species), ]

# Prepare the three matrices (one for each PC)
PC1_matrix <- as.matrix(PCs.plot[, "PC1", drop = FALSE])
PC2_matrix <- as.matrix(PCs.plot[, "PC2", drop = FALSE])
PC3_matrix <- as.matrix(PCs.plot[, "PC3", drop = FALSE])

rownames(PC1_matrix) <- rownames(PC2_matrix) <- rownames(PC3_matrix) <- PCs.plot$species

# List of matrices for the barplots
matrices <- list(PC1_matrix, PC2_matrix, PC3_matrix)

# Assign colors based on lifespan (0 = goldenrod1, 1 = darkorchid1)
lifespan_colors <- ifelse(PCs.plot$lifespan == 0, "goldenrod1", "darkorchid1")

# Set layout for plotting (1 row, 4 columns: 1 for tree + 3 for stacked barplots)
par(mfrow = c(1, 4))

# Step 1: Plot the tree and the first barplot (PC1)
plotTree.barplot(contree.rooted.PCA, matrices[[1]],
                 args.barplot = list(xlab = "PC1", col = lifespan_colors, mar = c(5.1, 0, 2.1, 2.1)), 
                 add = TRUE)

# Step 2: Add the second barplot (PC2) using add = TRUE, skip tree plotting
plotTree.barplot(contree.rooted.PCA, matrices[[2]],
                 args.barplot = list(xlab = "PC2", col = lifespan_colors, mar = c(5.1, 0, 2.1, 2.1)),
                 args.plotTree = list(plot = FALSE), add = TRUE)

# Step 3: Add the third barplot (PC3) using add = TRUE, skip tree plotting
plotTree.barplot(contree.rooted.PCA, matrices[[3]],
                 args.barplot = list(xlab = "PC3", col = lifespan_colors, mar = c(5.1, 0, 2.1, 2.1)),
                 args.plotTree = list(plot = FALSE), add = TRUE)

#########################
#### Stem resistance ####
#########################

supplement.BM1 <- mclapply(indices, function(i) {
  mvBM(simmaps.ls[[i]][[1]], 
       anatomy.scaled[,5:9], 
       model = "BM1", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(supplement.BM1, "./RDS/supplement.BM1.RDS")

supplement.BMM <- mclapply(indices, function(i) {
  mvBM(simmaps.ls[[i]][[1]], 
       anatomy.scaled[,5:9], 
       model = "BMM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(supplement.BMM, "./RDS/supplement.BMM.RDS")

supplement.OU1 <- mclapply(indices, function(i) {
  mvOU(simmaps.ls[[i]][[1]], 
       anatomy.scaled[,5:9], 
       model = "OU1", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(supplement.OU1, "./RDS/supplement.OU1.RDS")

supplement.OUM <- mclapply(indices, function(i) {
  mvOU(simmaps.ls[[i]][[1]], 
       anatomy.scaled[,5:9], 
       model = "OUM", 
       method = "rpf",
       scale.height = TRUE,
       control = list(maxit = 300000),
       param = list(root = FALSE))
}, mc.cores = 20)

saveRDS(supplement.OUM, "./RDS/supplement.OUM.RDS")

AICc.compare(supplement.BM1)
AICc.compare(supplement.BMM)
AICc.compare(supplement.OU1)
AICc.compare(supplement.OUM)

##############
#### AICc ####
##############

supplement.summary <- data.frame(Tree = indices,
                                 BM1 = AICc.vec(supplement.BM1),
                                 BMM = AICc.vec(supplement.BMM),
                                 OU1 = AICc.vec(supplement.OU1),
                                 OUM = AICc.vec(supplement.OUM))

supplement.summary.melt <- supplement.summary %>%
  gather(key = "Model", value = "AICc", -Tree) 

supplement.summary.melt <- supplement.summary.melt %>%
  group_by(Tree) %>%
  arrange(Tree, AICc)%>%
  mutate(ModelOrder = max(row_number()) + 1 - row_number(),
         Difference = case_when(
           row_number() == 1 ~ AICc,
           row_number() == 2 ~ AICc - lag(AICc),
           row_number() == 3 ~ AICc - lag(AICc, order_by = AICc, n = 1),
           row_number() == 4 ~ AICc - lag(AICc, order_by = AICc, n = 2)
         )
  )

supplement.summary.plot <- ggplot(supplement.summary.melt, aes(x = factor(Tree), y = Difference, fill = Model)) +
  geom_bar(stat = "identity", aes(group = ModelOrder)) +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_x_discrete(breaks = levels(factor(supplement.summary.melt$Tree))[seq(1, length(levels(factor(supplement.summary.melt$Tree))), 10)]) +
  scale_fill_manual(values = c("firebrick4", "indianred", "springgreen4", "palegreen")) +
  coord_cartesian(ylim = c(min(supplement.summary.melt$AICc), NA)) +
  theme(legend.position="none") +
  # scale_y_reverse() +
  labs(title = element_blank(),
       x = "Tree", y = "AICc", fill = "Model")

ggsave("./Plots/supplement_AICc.pdf", 
       supplement.summary.plot, 
       device = "pdf", 
       width = 6, 
       height = 1.5, 
       units = "in")

################
#### Thetas ####
################

supplement.theta <- theta.compare(supplement.OUM, anatoscale[5:9], anatomean[5:9])

supplement.theta.bind <- theta_converter(supplement.theta, "theta", "lifespan")

supplement.theta.bind <- supplement.theta.bind %>%
  mutate(column = factor(lifespan), 
         panel_id = factor(panel_id),
         list_id = factor(list_id))

supplement.theta.bind$lifespan <- gsub("0", "annuals", supplement.theta.bind$lifespan)
supplement.theta.bind$lifespan <- gsub("1", "perennials", supplement.theta.bind$lifespan)

supplement.theta.plot <- ggplot(supplement.theta.bind, aes(x = value, fill = interaction(list_id, lifespan))) + 
  geom_density(alpha = 0.6) + # Add some transparency
  facet_wrap(~row_name, scales = "free", ncol = 1) + # Create a panel for each row, adjust ncol as needed
  labs(title = element_blank(), 
       x = "", 
       y = "Density") +
  theme_minimal() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_fill_manual(values = c("goldenrod1", "darkorchid1")) +
  theme(legend.position = "none") +
  labs(fill = "")

ggsave("./Plots/supplement_OUM_theta.pdf", 
       supplement.theta.plot, 
       device = "pdf", 
       width = 6, 
       height = 6, 
       units = "in")