#Function returning mean and sd from AICc scores

AICc.compare <- function (data) {
  AICc <- list()
  for (i in seq_along(data)){
    AICc[[i]] <- data[[i]]$AICc
  }
  return (c (mean(unlist(AICc)), sd(unlist(AICc))))
} 

AICc.vec <- function (data) {
  AICcvec <- list()
  for (i in seq_along(data)){
    AICcvec[[i]] <- data[[i]]$AICc
  }
  return (unlist(AICcvec))
} 

GICvec <- function (data) {
  GIC <- list()
  for (i in seq_along(data)){
    GIC[[i]] <- data[[i]]$GIC
  }
  return (unlist(GIC))
} 

EICvec <- function (data) {
  EIC <- list()
  for (i in seq_along(data)){
    EIC[[i]] <- data[[i]][[1]]
  }
  return (unlist(EIC))
} 

pvalvec <- function (data, x) {
  pval <- list()
  for (i in seq_along(data)){
    pval[[i]] <- data[[i]]$pvalue[[x]]
  }
  return (unlist(pval))
} 

theta.compare <- function (data, scale, center) {
  theta <- list()
  for (i in seq_along(data)){
    theta[[i]] <- apply (data[[i]]$theta, 1, function(x) x * scale + center)
  }
  return (theta)
} 

theta.compare.BIO <- function (data) {
  theta <- list()
  for (i in seq_along(data)){
    theta[[i]] <- t(data[[i]]$theta)
  }
  return (theta)
} 

theta_converter <- function(df_list, list_id, names) {
  bind_rows(lapply(seq_along(df_list), function(i) {
    as.data.frame(df_list[[i]]) %>%
      mutate(row_name = rownames(.),
             df_id = i,
             list_id = list_id) %>%
      pivot_longer(cols = c('0', '1'), names_to = names, values_to = "value")
  }), .id = "panel_id")
}

theta_converter_0 <- function(df_list, list_id, names) {
  bind_rows(lapply(seq_along(df_list), function(i) {
    as.data.frame(df_list[[i]]) %>%
      mutate(row_name = rownames(.),
             df_id = i,
             list_id = list_id) %>%
      pivot_longer(cols = c('OU1'), names_to = names, values_to = "value")
  }), .id = "panel_id")
}

pvalvec_phylolm <- function (data, j) {
  pval <- list()
  for (i in seq_along(data)){
    pval[[i]] <- summary(data[[i]])[[2]][j,6]
  }
  return (unlist(pval))
}

#Function calculating mean Q matrix

Q.avg <- function (corHMM) {
  Q <- list()
  for (i in seq_along(corHMM)){
    Q[[i]] <- corHMM[[i]]$solution
  }
  Q.dcb <- do.call (cbind, Q)
  Q.array <- array (Q.dcb, 
                    dim = c(dim(Q[[1]]), length(Q)))
  Q.mean <- apply (Q.array, 
                   c(1,2), 
                   mean, 
                   na.rm = TRUE)
  return(Q)
  return(Q.mean)
} 

#Modifiet corHMM plot function
plotRECONJB <- function (phy, likelihoods, tipcolors, piecolors = NULL, cex = 0.5, pie.cex = 0.25, 
                         file = NULL, height = 11, width = 8.5, show.tip.label = TRUE, 
                         title = NULL, ...) 
{
  if (is.null(piecolors)) {
    piecolors = c("white", "black", "red", "yellow", "forestgreen", 
                  "blue", "coral", "aquamarine", "darkorchid", "gold", 
                  "grey", "yellow", "#3288BD", "#E31A1C")
  }
  if (!is.null(file)) {
    pdf(file, height = height, width = width, useDingbats = FALSE)
  }
  plot(ladderize(phy), cex = cex, show.tip.label = show.tip.label, tip.color = tipcolors, ...)
  if (!is.null(title)) {
    title(main = title)
  }
  nodelabels(pie = likelihoods, piecol = piecolors, cex = pie.cex)
  states <- colnames(likelihoods)
  if (!is.null(file)) {
    dev.off()
  }
}

#Summarizing model fitting
summarize_model_fit <- function(data) {
  best_models <- data %>% filter(ModelOrder == 4)
  second_best_models <- data %>% filter(ModelOrder == 3)
  
  combined <- best_models %>% 
    inner_join(second_best_models, by = "Tree", suffix = c("_best", "_second_best"))
  
  combined <- combined %>%
    mutate(Difference = AICc_second_best - AICc_best)
  
  best_fit_counts <- combined %>%
    mutate(FitType = ifelse(Difference > 2, "BestFit", "Inconclusive")) %>%
    group_by(Model_best, FitType) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = FitType, values_from = Count, values_fill = list(Count = 0))
  
  summary_table <- best_fit_counts %>%
    mutate(Total = BestFit + Inconclusive)

}

#Wrap splitter
save_facet_plots <- function(plot, facet_var, data, path_prefix) {
  facets <- unique(data[[facet_var]])
  
  for (facet in facets) {
    p_facet <- plot +
      ggforce::facet_wrap_paginate(~row_name, scales = "free", ncol = 1, page = which(facets == facet)) +
      labs(title = paste("Hypothesis 4:", facet))
    
    ggsave(filename = paste0(path_prefix, "_", facet, ".tiff"), plot = p_facet, device = "tiff")
  }
}

