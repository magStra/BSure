#' Extracts core output from the output file created by the BSure function.
#' The extracted output is used in postprocessing functions for plotting and quality control.
#' @param BSure_output  output created by BSure function
#' @title extract_from_output
#' @return sample_summary with mean estimates (mean) and 95% credible intervals (CI_lower_bounds
#' for the lower bounds of the intervals, CI_upper bounds for the upper ones) of the essentiality
#' level of each gene, the probability of a gene being highly essential (prob_essential_II; the posterior distribution
#' of the essentiality being shifted less than 1/3 compared to a typical core-essential gene), the
#' probability of a gene being essential of level I (no overlap with distribution of typical nonessential gene),
#' and convergence statistics (Rhat=Gelman-Rubin R-hat statistic, tail_ESS=tail essential sample size)
#' In addition to the sample summary: prop_expensive_sampling: proportion of genes which
#' required a more expensive sampling algorithm
#' prop_very_expensive_sampling: proportion of genes which required the most expensive sampling
#' algorithm
#' @export
extract_from_output <- function(BSure_output)
{
  ncl <- ncol(BSure_output)
  res <- c()
  Rhat <- rep(0,ncl)
  CI_lower_bounds <- rep(0,ncl)
  CI_upper_bounds <- rep(0,ncl)
  prob_essential_I <- rep(0,ncl)
  prob_essential_II <- rep(0,ncl)
  prob_NE_05 <- rep(0,ncl)
  prob_essential_05 <- rep(0,ncl)
  tail_ESS <- rep(0,ncl)
  means <- rep(0,ncl)
  gene_names <- c()
  nExpensiveSampling <- 0
  nVeryExpensiveSampling <- 0
  for (j in 1:ncl)
  {
    Rhat[j] <- BSure_output[,j]$rhat[1]
    nExpensiveSampling = nExpensiveSampling+as.numeric(BSure_output[,j]$expensiveSampling)
    nVeryExpensiveSampling = nVeryExpensiveSampling+as.numeric(BSure_output[,j]$veryExpensiveSampling)
    CI_lower_bounds[j] <- BSure_output[,j]$quant025[1]
    CI_upper_bounds[j] <- BSure_output[,j]$quant975[1]
    tail_ESS[j] <- BSure_output[,j]$ess_tail[1]
    means[j] <- BSure_output[,j]$mean[1]
    prob_essential_II[j] <- BSure_output[,j]$probability_essential_II
    prob_essential_I[j] <- BSure_output[,j]$probability_essential_I
    prob_NE_05 <- BSure_output[,j]$probability_NE_05
    prob_essential_05 <-  BSure_output[,j]$probability_essential_05
    gene_names <- c(gene_names,BSure_output[,j]$gene_names)
  }
  prop_expensive_sampling <- nExpensiveSampling/ncl
  prop_very_expensive_sampling <- nVeryExpensiveSampling/ncl
  sample_summary <- data.frame(means=means,CI_lower_bounds=CI_lower_bounds,CI_upper_bounds=CI_upper_bounds,
                              Rhat=Rhat,tail_ESS=tail_ESS,prob_essential_II=prob_essential_II,
                              prob_essential_I=prob_essential_I,
                              prob_NE_05=prob_NE_05,prob_essential_05=prob_essential_05,gene_names=gene_names)
  res$prop_expensive_sampling <- prop_expensive_sampling
  res$prop_very_expensive_sampling <- prop_very_expensive_sampling
  res$sample_summary <- sample_summary
  return(res)
}

#' This function extracts sets of
#' 1) highly essential genes (essentiality level II, the posterior distribution
#' of the essentiality is shifted less than 1/3 compared to a typical core-essential gene
#' with with FDR of 0.05, unless a different FDR is specified
#' 2) genes confidently not nonessential with FDR of 0.01, unless a different FDR is specified (essentiality level I)
#' @title extract_gene_categories
#' @import eulerr
#' @param extracted_output output of the function extract_from_output
#' @param FDR_II set to 0.05 by default; false positive rate for essentiality of type II
#' @param FDR_I false positive rate for essentiality of type I; set to 0.05 by default
#' @return essential_genes_II  genes of essentiality level II (see above)
#' essential_genes_I genes of essentiality level I (see above).
#' proportion_core_essential_in_essential_II which proportion of core essential genes is in the essential genes
#' II group
#' (as core essential genes we use the intersection of core essential genes from  Hart and Moffat (2016) and
#' Behan et al. (2019))
#' proportion_core_essential_in_essential_I proportion of core essential genes in essential I group
#' proportion_non_essential_genes_in_essential_II proportion of nonessential genes in essential II group
#' (nonessential as in Hart and Moffat (2016)); this number should be zero or very close to 0
#' proportion_non_essential_genes_in_essential_I proportion of nonessential genes in essential I group
#' score A summary score for the ability of a screen to separate core essential and nonessential genes.
#'
#' @references Hart T, Moffat J.
#' BAGEL: a computational framework for identifying essential genesfrom pooled library screens.
#' BMC Bioinformatics. 2016;17(1):164.  doi:10.1186/s12859-016-1015-8
#' Behan FM et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens.
#' Nature. 2019;568(7753):511–516
#' @export
extract_gene_categories <- function(extracted_output,FDR_II=0.05,FDR_I= 0.05,figures=T)
{
  data(auxData)
  coreEssGenes <- intersect(CEGv2,coreFitnessBehanEtAl)
  nonEssGenes <- NEG
  output <- c()
  coreEssGenes <- intersect(coreEssGenes,extracted_output$sample_summary$gene_names)
  nonEssGenes <- intersect(nonEssGenes,extracted_output$sample_summary$gene_names)
  probs_essential_II <- extracted_output$sample_summary$prob_essential_II
  alpha0 <- .findOptimalAlpha(probs_essential_II,FDR_II)
  probs_essential_I <- extracted_output$sample_summary$prob_essential_I
  alpha1 <-.findOptimalAlpha(probs_essential_I,FDR_I)
  output$essential_genes_II<- extracted_output$sample_summary$gene_names[probs_essential_II>alpha0]
  output$essential_genes_I <- extracted_output$sample_summary$gene_names[probs_essential_I>alpha1]
  coreEssGenesDiscovered <- intersect(output$essential_genes_II,coreEssGenes)
  coreEssGenesDiscoveredNotNE <- intersect(output$essential_genes_I,coreEssGenes)
  output$proportion_core_essential_in_essential_II <- length(coreEssGenesDiscovered)/length(coreEssGenes)
  output$proportion_core_essential_in_essential_I <- length(coreEssGenesDiscoveredNotNE)/length(coreEssGenes)
  non_essential_genes_in_essential_II <- intersect(output$essential_genes_II,nonEssGenes)
  non_essential_genes_in_essential_I <- intersect(output$essential_genes_I,nonEssGenes)
  output$proportion_non_essential_genes_in_essential_II <- length(non_essential_genes_in_essential_II)/length(nonEssGenes)
  output$proportion_non_essential_genes_in_essential_I <- length(non_essential_genes_in_essential_I)/length(nonEssGenes)
  gene_set_list <- list(output$essential_genes_II,coreEssGenes,nonEssGenes)
  names_set <- c("essential (II)",
                 "core essential genes",
                 "nonessential genes")
  names(gene_set_list) <- names_set
  if(figures==T)
  print(.plot_Venn_diagram(gene_set_list))
  gene_set_list <- list(output$essential_genes_I,coreEssGenes,nonEssGenes)
  names_set <- c("essential (I)",
                 "core essential genes",
                 "nonessential genes")
  names(gene_set_list) <- names_set
  if(figures==T)
  print(.plot_Venn_diagram(gene_set_list))
  gene_set_list <- list(output$essential_genes_I,output$essential_genes_II)
  names_set <- c("essential (I)","essential (II)")
  names(gene_set_list) <- names_set
  if(figures==T)
  print(.plot_Venn_diagram(gene_set_list))
  output$score <- sqrt(output$proportion_core_essential_in_essential_II*output$proportion_core_essential_in_essential_I)-
    output$proportion_non_essential_genes_in_essential_I
  return(output)
}

#' This function is identical to extract_gene_categories except it finds sets of essential genes given a false
#' negative rather than a false positive rate.
#' @title extract_gene_categories_FNR
#' @param extracted_output output of the function extract_from_output
#' @param FNR_II set to 0.05 by default; false negative rate for essentiality of type II
#' @param FNR_I set to 0.05 by default
#' @return essential_genes_II  genes of essentiality level II
#' essential_genes_I genes of essentiality level I
#' @references Hart T, Moffat J.
#' BAGEL: a computational framework for identifying essential genesfrom pooled library screens.
#' BMC Bioinformatics. 2016;17(1):164.  doi:10.1186/s12859-016-1015-8
#' Behan FM et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens.
#' Nature. 2019;568(7753):511–516
#' @export
extract_gene_categories_FNR <- function(extracted_output,FNR_II=0.05,FNR_I= 0.05)
{
  data(auxData)
  coreEssGenes <- intersect(CEGv2,coreFitnessBehanEtAl)
  nonEssGenes <- NEG
  output <- c()
  coreEssGenes <- intersect(coreEssGenes,extracted_output$sample_summary$gene_names)
  nonEssGenes <- intersect(nonEssGenes,extracted_output$sample_summary$gene_names)
  probs_essential_II <- extracted_output$sample_summary$prob_essential_II
  alpha0 <- .findOptimalAlpha_FNR(probs_essential_II,FNR_II)
  probs_essential_I <- extracted_output$sample_summary$prob_essential_I
  alpha1 <-.findOptimalAlpha_FNR(probs_essential_I,FNR_I)
  output$essential_genes_II<- extracted_output$sample_summary$gene_names[probs_essential_II>alpha0]
  output$essential_genes_I <- extracted_output$sample_summary$gene_names[probs_essential_I>alpha1]
  return(output)
}


library(RColorBrewer)

.false_positive_rate <- function(probs,alpha)
{
  return(sum((1-probs)*(probs>alpha))/sum(probs>alpha))
}

.false_negative_rate <- function(probs,alpha)
{
  return(sum(probs*(probs<=alpha))/sum(probs<=alpha))
}


#' @import rootSolve
.findOptimalAlpha <- function(probs,FDR=0.05)
{
  if (.false_positive_rate(probs,max(probs)-0.0001) >= FDR)
    return(max(probs)-0.0001) else{
    return(uniroot(function(alpha) return(.false_positive_rate(probs,alpha)-FDR),c(min(probs)+0.0001,max(probs)-0.0001))$root)}}

.findOptimalAlpha_FNR <- function(probs,FNR=0.05)
{ if (.false_negative_rate(probs,min(probs)+0.0001) >= FNR)
  return(min(probs)+0.0001)
  else if (.false_negative_rate(probs,max(probs)-0.0001) <= FNR) {return(max(probs)-0.0001)}else{
  return(min(1,uniroot(function(alpha) return(.false_negative_rate(probs,alpha)-FNR),c(min(probs)+0.0001,max(probs)-0.0001))$root))}
}

#' @import eulerr
.plot_Venn_diagram <- function(gene_set_list,file_name,res=1200)
{
  p <- plot(euler(gene_set_list, shape = "ellipse"), quantities = TRUE,adjust_labels=TRUE)
  return(p)

}

#' @title find_differential_genes
#' Compares to screens and finds genes that are essential of type I or II in the first screen (controlling for FDR=0.001),
#' but not in the second (controlling for both FNR=FNR_II and FDR=FDR_II, to ensure low rate of false positives of differential essentiality)
#' This function may only be used if the second screen has sufficient power to detect essential genes compared to the first screen.
#' If this balance of power is lacking, calling the function will results in an error.
#' @params extracted_output_1 output of the function extract_from_output for the first screen
#' @params extracted_output_2 output of the function extract_from_output for the second screen
#' @params FDR_II FDR threshold for both first and second screen
#' @params FNR_II FNR threshold for second screen
#' @return  genes that are essential of type II in the first screen (controlling for FDR),
#' but not in the second (controlling for FNR)
#' @export


find_differential_genes <- function(extracted_output_1,extracted_output_2,
                                    FDR_II=0.05,FNR_II=0.1,FNR_I=0.05)
{
  A <- compare_to_gold_standard(extracted_output_1,extracted_output_2)$proportion_genes_recovered_II_FNR
  if (A < 0.75)
  {
    print("Less than 75% of genes essential of type II in screen one at the given FDR are essential of type II in screen 2 at the given FNR. This indicates a power imbalance across the screens to detect essential genes, and therefore this method to detect differential essentiality cannot be used.")
    return(NULL)
  }
  extracted_genes_1 <- extract_gene_categories(extracted_output_1,FDR_II=FDR_II,figures=F)
  extracted_genes_2 <- extract_gene_categories_FNR(extracted_output_2,FNR_II=FNR_II)
  extracted_genes_3 <- extract_gene_categories(extracted_output_2,FDR_II=FDR_II,figures=F)
  output <- c()
  genes_1_II <- as.character(extracted_genes_1$essential_genes_II)
  genes_2 <- unique(c(as.character(extracted_genes_2$essential_genes_II),
                      as.character(extracted_genes_3$essential_genes_II)))
  output <- setdiff(genes_1_II,genes_2)
  return(output)
}

#' @title compare_to_gold_standard
#' Compares a new screen to your gold-standard screen. For example, your gold-standard screen may be of high
#' coverage and your new screen of lower coverage, and you want to check whether essential genes of type I or II
#' from the gold standard screen are also recovered by the new screen.
#' @params extracted_output_gold_standard output of the function extract_from_output for the gold standard screen
#' @params extracted_output_new output of the function extract_from_output for the new screen
#' @return proportion_genes_recovered_I Proportion of the genes essential of type I in the first screen that are
#' also essential of type I in the second screen.
#' @return proportion_genes_recovered_II Proportion of the genes essential of type II in the first screen that are
#' also essential of type I in the second screen.
#' @return proportion_genes_recovered_I_FNR Proportion of the genes essential of type I in the first screen that are
#' also essential of type I in the second screen, when FNR is controlled at 0.05 in the second screen. To avoid cases with
#' very low FNR we use the union of the FNR and FDR controlled sets.
#' @return proportion_genes_recovered_II_FNR Proportion of the genes essential of type II in the first screen that are
#' also essential of type I in the second screen, when FNR is controlled at 0.05 in the second screen. Again we use
#' the union of the FNR and FDR controlled sets.
#' @export

compare_to_gold_standard <- function(extracted_output_gold_standard,extracted_output_new)
{
  extracted_genes_1 <- extract_gene_categories(extracted_output_gold_standard,figures=F,FDR_I=0.05,FDR_II=0.05)
  extracted_genes_2 <- extract_gene_categories_FNR(extracted_output_new,FNR_I=0.05,FNR_II=0.05)
  extracted_genes_3 <- extract_gene_categories(extracted_output_new,figures=F,FDR_I=0.05,FDR_II=0.05)
  output <- c()
  genes_1_I <- as.character(extracted_genes_1$essential_genes_I)
  genes_1_II <- as.character(extracted_genes_1$essential_genes_II)
  genes_2a <- as.character(extracted_genes_3$essential_genes_II)
  genes_3a <- as.character(extracted_genes_3$essential_genes_I)
  genes_2 <- unique(c(as.character(extracted_genes_2$essential_genes_II),genes_2a))
  genes_3 <- unique(c(as.character(extracted_genes_2$essential_genes_I),genes_3a))
  output$proportion_genes_recovered_II_FNR <- length(intersect(genes_1_II,genes_2))/length(genes_1_II)
  output$proportion_genes_recovered_I_FNR <- length(intersect(genes_1_I,genes_3))/length(genes_1_I)
  output$proportion_genes_recovered_II <- length(intersect(genes_1_II,genes_2a))/length(genes_1_II)
  output$proportion_genes_recovered_I <- length(intersect(genes_1_I,genes_3a))/length(genes_1_I)
  return(output)
}


#' @title compare_length_of_credible_intervals
#' Compares the length of the credible intervals of the gene essentiality estimates for two different screens or
#' studies.
#' @param extracted_output_1 extracted output for screen/study 1, computed from BSure output using the function
#' extract_from_output.
#' @param extracted_output_2 extracted output for screen/study 2, computed from BSure output using the function
#' extract_from_output.
#' @param study_name_1 name of first screen/study
#' @param study_name_2 name of second screen/study
#' @return ratio of credible intervals, prints scatterplot comparing lengths of credible intervals.
#' @export

compare_length_of_credible_intervals <- function(extracted_output_1,extracted_output_2,study_name_1,study_name_2)
{
  CIs_1 <- extracted_output_1$sample_summary$CI_upper_bounds -
    extracted_output_1$sample_summary$CI_lower_bounds
  CIs_2 <- extracted_output_2$sample_summary$CI_upper_bounds -
    extracted_output_2$sample_summary$CI_lower_bounds
  df <- data.frame(CIs_1=CIs_1,CIs_2=CIs_2)
  pp <- ggplot(df,aes(y=CIs_1,x=CIs_2)) + geom_point(alpha=0.05,size=1,color="blue")+
    ylab(study_name_1)+xlab(study_name_2)+theme_classic()+theme(axis.text.x = element_text(size=15),
                                                            axis.text.y = element_text(size=15),axis.title.x = element_text(size=20),
                                                            axis.title.y=element_text(size=20),plot.title = element_text(size=20))+
    ggtitle("length of credible interval")+ geom_abline(aes(intercept = 0, slope = 1),
                                                     color="black",linetype = "dashed",size=1.2,show.legend = T)
  print(pp,width=8.8/2.54)
  output <- cbind(CIs_1,CIs_2/CIs_1)
  colnames(output) <- c("study 1","study 2/study 1")
  return(cbind(CIs_1,CIs_2/CIs_1))
}
#'
#'
#'
#'
#'
#'
#' @title plot_credible_intervals
#' @import ggplot2
#' @import ggthemes
#' @param extracted_output output of the function extract_from_output
#' @param name name of the screen for naming of output folder and pdf files
#' @return creates pdf files with the following plots:
#' 1) credible_intervals_as_function_of_mean: credible intervals of essentiality estimates along y-axis,
#' mean of essentality estimates on x-axis
#' 2) credible_intervals_ordered_by_mean: credible intervals of essentiality estimates along y-axis,
#' ordered by mean of essentiality estimates (ordered along x-axis)
#' 3) density_length_of_credible_interval: density distribution of the lengths of the credible intervals
#' 4) tail_ESS_credible_interval and Rhat_credible_interval: diagnostic plots to identify any correlation
#' between convergence measures and the lengths of credible intervals. If there is correlation (which is rare)
#' this may point to convergence problems of the sampler, and the results may be less reliable.
#' @export
plot_credible_intervals <- function(extracted_output,name)
{
  dir.create(paste0(name,"/"))
  name <- paste0(name,"/",name)
  mu <- extracted_output$sample_summary$means
  names(mu) <- extracted_output$gene_names
  minCredInt <- extracted_output$sample_summary$CI_lower_bounds
  names(minCredInt) <- extracted_output$gene_names
  maxCredInt <- extracted_output$sample_summary$CI_upper_bounds
  names(maxCredInt) <- extracted_output$gene_names
  tail_ESS <-  extracted_output$sample_summary$tail_ESS
  names(tail_ESS) <- extracted_output$gene_names
  Rhat <- extracted_output$sample_summary$Rhat
  names(Rhat) <- extracted_output$gene_names
  aa <- order(mu)
  boundCredInt <- c(minCredInt[aa],maxCredInt[aa])
  postMean <- c(mu[aa],mu[aa])
  minMaxInd <- rep(c("lower","upper bound"),each=length(mu))
  df <- data.frame(postMean=postMean,boundCredInt = boundCredInt, minMaxInd=minMaxInd)
  pp<-ggplot(df, aes(x=postMean,y=boundCredInt)) +
    geom_line(aes(color = minMaxInd)) +
    theme_bw(base_size = 11)+xlab("posterior mean")+
    ylab("95% credible intervals")+scale_color_colorblind()+theme(legend.title = element_blank(),
                                                               legend.position = "top",legend.text = element_text(size=13),legend.title.align = 0)
  pdf(file=paste(name,"_credible_intervals_as_function_of_mean.pdf",sep=""),width=unit(3.5,"cm"),height=unit(2.5,"cm"))
  print(pp,dpi=1200)
  dev.off()
  ordering = rep(1:length(aa),2)
  df <- data.frame(ordering=ordering,boundCredInt = boundCredInt, minMaxInd=minMaxInd)
  pp<-ggplot(df, aes(x=ordering,y=boundCredInt)) +
    geom_line(aes(color = minMaxInd)) +
    theme_bw(base_size = 11)+xlab("rank of mean essentiality")+
    ylab("95% credible intervals")+scale_color_colorblind()+theme(legend.title = element_blank(),
                                                               legend.position = "top",legend.text = element_text(size=13),legend.title.align = 0)
  pdf(file=paste(name,"_credible_intervals_ordered_by_mean.pdf",sep=""),width=unit(3.5,"cm"),height=unit(2.5,"cm"))
  print(pp,dpi=1200)
  dev.off()
  lengthCredInt <- maxCredInt-minCredInt
  df <- data.frame(lengthCredInt=lengthCredInt)
  pp<-ggplot(df,aes(x=lengthCredInt)) + geom_density(position=position_stack(vjust=0.5))+theme_classic(base_size = 11)+xlab("Length of credible interval")+
    ylab("number of genes")+geom_rug(color="blue",alpha=0.1)+theme(axis.text.x = element_text(size=11),
             axis.text.y = element_text(size=11))
  pdf(file=paste(name,"_density_length_of_credible_interval.pdf",sep=""),width=unit(4,"cm"),height=unit(2.5,"cm"))
  print(pp,dpi=1200)
  dev.off()
  df <- data.frame(tail_ESS=tail_ESS,lengthCredInt=lengthCredInt)
  pp <- ggplot(df,aes(x=tail_ESS,y=lengthCredInt))+geom_point()+xlab("tail ESS")+ylab("length of 95%\n credible interval")+theme_classic(base_size = 11)+
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11))
  pdf(file=paste(name,"_tail_ESS_credible_interval.pdf",sep=""),width=unit(4,"cm"),height=unit(2.5,"cm"))
  print(pp,dpi=1200)
  dev.off()
  df <- data.frame(Rhat=Rhat,lengthCredInt=lengthCredInt)
  pp <- ggplot(df,aes(x=Rhat,y=lengthCredInt))+geom_point()+xlab("R-hat")+ylab("length of 95%\n credible interval")+theme_classic(base_size = 11)+
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11))
  pdf(file=paste(name,"_Rhat_credible_interval.pdf",sep=""),width=unit(4,"cm"),height=unit(2.5,"cm"))
  print(pp,dpi=1200)
  dev.off()
}
#'
#'
#' technical_replicate_correlations_essential <- function(lfc,plotTitle="")
#' {
#'   l <- floor(nrow(lfc)*0.1)
#'   xx <- complete.cases(lfc)
#'   A <- lfc[xx,-(1:2)]
#'   x1 <- order(rowMeans(A),decreasing = F)[1:l]
#'   corMatrix <- cor(A[x1,])
#'   melted_corMatrix <- melt(corMatrix)
#'   pRepCor <- ggplot(data = melted_corMatrix, aes(x=Var1, y=Var2, fill=value)) +
#'     geom_tile()+theme_minimal()+scale_fill_gradient(name="Correlation\n of log-fold changes\n 20% most depleted")+
#'     theme(legend.key.height = unit(0.4, "cm"),legend.text=element_text(size=14),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=12),
#'           legend.title = element_text(size=18),axis.title = element_text(size=14),plot.title=element_text(size=14),axis.text.y = element_text(size=12)) +
#'     ylab("")+xlab("")+ ggtitle(plotTitle)
#'   output <- c()
#'   output$correlation_matrix <- corMatrix
#'   output$plot <- pRepCor
#'   output$meanCorrelation <- mean(corMatrix[upper.tri(corMatrix)])
#'   return(output)
#' }
#'
#' technical_replicate_correlations <- function(lfc,plotTitle="")
#' {
#'   #compute a score so that we can compare this across screens
#'   #and make plots
#'   lfc1 <- lfc[complete.cases(lfc[,-(1:2)]),-c(1:2)]
#'   corMatrix <- cor(lfc1)
#'   melted_corMatrix <- melt(corMatrix)
#'   pRepCor <- ggplot(data = melted_corMatrix, aes(x=Var1, y=Var2, fill=value)) +
#'     geom_tile()+theme_minimal()+scale_fill_gradient(name="Correlation\n of log-fold changes")+
#'     theme(legend.key.height = unit(0.4, "cm"),legend.text=element_text(size=14),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size=12),
#'           legend.title = element_text(size=18),axis.title = element_text(size=14),plot.title=element_text(size=14),axis.text.y = element_text(size=12)) +
#'     ylab("")+xlab("")+ ggtitle(plotTitle)
#'   output <- c()
#'   output$correlation_matrix <- corMatrix
#'   output$plot <- pRepCor
#'   output$meanCorrelation <- mean(corMatrix[upper.tri(corMatrix)])
#'   return(output)
#' }
#'
#' biological_replicate_correlations <- function(lfcRep1,lfcRep2)
#' {
#'   xx <- complete.cases(lfcRep1)
#'   xx[!(complete.cases(lfcRep2))] <- F
#'   return(mean(cor(lfcRep1[xx,-(1:2)],lfcRep2[xx,-(1:2)])))
#' }
#'
#' biological_replicate_correlations_essential <- function(lfcRep1,lfcRep2)
#' {
#'   l <- floor(nrow(lfcRep1)*0.2)
#'   xx <- complete.cases(lfcRep1)
#'   xx[!(complete.cases(lfcRep2))] <- F
#'   A <- lfcRep1[xx,-(1:2)]
#'   B <- lfcRep2[xx,-(1:2)]
#'   x1 <- order(rowMeans(cbind(A,B)),decreasing = F)[1:l]
#'
#'   cor(A[x1,],B[x1,])
#'   mean(cor(A[x1,],B[x1,]))
#'   return(mean(cor(A[x1,],B[x1,])))
#' }
#'
#'
#'
