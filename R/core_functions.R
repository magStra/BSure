#' Computes log-fold changes of gRNA counts from gRNA count and control data.
#' @title compute_lfc_from_counts
#' @param counts MxR matrix of gRNA counts, where M is the number of gRNAs and R the number of replicates
#' @param controls Mx1 vector or MxC matrix, if there are C different control measurements
#' @return matrix of log-fold changes
#' @author Magdalena Strauss
#' @example data(HT29_count_data)
#' counts <- HT29_count_data[,-c(1,2,ncol(HT29_count_data))]
#' controls <- HT29_count_data$HT29_controls
#' lfc_HT29 <- compute_lfc_from_counts(counts,controls)
#' @export
compute_lfc_from_counts <- function(counts,controls)
{
  colSumsCounts <- colMeans(counts,na.rm = T)*nrow(counts)
  controls <- as.matrix(controls)
  counts <- as.matrix(counts)
  counts <- t(apply(counts,1,function(x) {return(x/colSumsCounts*10^6)}))
  colSumsControls <- colMeans(controls,na.rm = T)*nrow(controls)
  controls <- t(apply(controls,1,function(x) {return(x/colSumsControls*10^6)}))
  v <- rowMeans(log2(controls+1),na.rm=T)
  diff <- apply(as.matrix(counts),2,function(x) return(log2(x+1)-v))
  rownames(diff) <- rownames(counts)
  return((diff))
}


#' @title compute_mean_lfc
#' @param lfc_matrix MxR matrix of log-fold changes, where M is the number of gRNAs and R the number of replicates
#' @param gene_vector vector of length M listing the genes corresponding to the gRNAs (in the same order as the rows of lfc)
#' @return vector of mean log-fold changes across replicates
#' @export
#' @author Magdalena Strauss

compute_mean_lfc <- function(lfc_matrix,gene_vector)
{
  uniqueGeneVector <- unique(gene_vector)
  nGenes <- length(uniqueGeneVector)
  meanLfc <- rep(0,nGenes)
  sdsLfc <- rep(0,nGenes)
  for (j in 1:nGenes)
  {
    meanLfc[j] <- mean(as.vector(lfc_matrix[gene_vector == uniqueGeneVector[j],]),na.rm=T)
    sdsLfc[j] <- sd(as.vector(lfc_matrix[gene_vector == uniqueGeneVector[j],]),na.rm=T)
  }
  output <- c()
  output$mean_lfc <- meanLfc
  output$sds_lfc <- sdsLfc
  output$genes <- uniqueGeneVector
  return(output)
}



#' @title fit_bivariate_nv2
#' @param mean_essential vector of mean log-fold changes for all gRNAs targeting core essential genes,
#' where the mean is taken across replicates or cell lines and replicates, if there are several cell lines
#' @param mean_NE vector of mean log-fold changes for all gRNAs targeting listed nonessential genes
#' @param sd_essential standard deviation of log-fold changes for all gRNAs targeting core essential genes
#' @param sd_NE standard deviation of log-fold changes for all gRNAs targeting listed nonessential genes
#' @param plot_file name (without extension) of the file with the plot of the data and prior distribution
#' @return parameters of bivariate normal prior distribution for mean and standard deviation of gene essentiality
#' @export
#' @author Magdalena Strauss

#' @import ggplot2
#' @import ggpubr
#' @import MASS
#' @importFrom rstan optimizing
fit_bivariate_nv2 <- function(mean_essential,mean_NE,sd_essential,sd_NE,plot_file_name=NA)
{
  output <- c()
  xx <- !(is.na(sd_essential))
  mean_essential <- mean_essential[xx]
  sd_essential <- sd_essential[xx]
  xx <- !(is.na(sd_NE))
  mean_NE <- mean_NE[xx]
  sd_NE <- sd_NE[xx]
  output$essential <- optimizing(stanmodels$computeParamsMultivariateN,data=list(M=length(mean_essential),y=cbind(mean_essential,sd_essential)))
  output$NE <- optimizing(stanmodels$computeParamsMultivariateN,data=list(M=length(mean_NE),y=cbind(mean_NE,sd_NE)))
  if (!(is.na(plot_file_name)))
  {dfE1 <- data.frame(mean_essential,sd_essential)
  pE1 <- ggplot(dfE1, aes(x=mean_essential,y=sd_essential) ) +theme(axis.text.x = element_text(size=13),
                                                                    axis.text.y = element_text(size=13),axis.title.x = element_text(size=12),
                                                                    axis.title.y=element_text(size=12))+
    geom_bin2d(bins=10)+ xlab("mean") + ylab("sd")+
    scale_fill_continuous(type = "viridis") + theme_bw()+ ggtitle("Core essential - data")
  dfNE1 <- data.frame(mean_NE,sd_NE)
  pNE1 <- ggplot(dfNE1, aes(x=mean_NE,y=sd_NE) ) +theme(axis.text.x = element_text(size=13),
                                                        axis.text.y = element_text(size=13),axis.title.x = element_text(size=12),
                                                        axis.title.y=element_text(size=12))+
    geom_bin2d(bins=10) + xlab("mean") + ylab("sd")+
    scale_fill_continuous(type = "viridis") + theme_bw()+ ggtitle("Nonessential - data")}
  SigmaEssential <- matrix(output$essential$par[9:12],nrow=2)
  SigmaNE <- matrix(output$NE$par[9:12],nrow=2)
  mean_essential <- output$essential$par[5:6]
  mean_NE <- output$NE$par[5:6]
  samplesEsssential <- mvrnorm(n = 500, mean_essential, SigmaEssential)
  samplesmean_essential <- samplesEsssential[,1]
  samplessd_essential <- samplesEsssential[,2]
  samplesNE <- mvrnorm(n = 500, mean_NE, SigmaNE)
  samplesmean_NE <- samplesNE[,1]
  samplessd_NE <- samplesNE[,2]
  if (!(is.na(plot_file_name)))
  {dfE2 <- data.frame(samplesmean_essential=samplesmean_essential,samplessd_essential=samplessd_essential)
  dfNE2 <- data.frame(samplesmean_NE=samplesmean_NE,samplessd_NE=samplessd_NE)
  pE2 <- ggplot(dfE2, aes(x=samplesmean_essential,y=samplessd_essential) ) +theme(axis.text.x = element_text(size=13),
                                                                                  axis.text.y = element_text(size=13),axis.title.x = element_text(size=12),
                                                                                  axis.title.y=element_text(size=12))+
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) + xlab("mean") + ylab("sd")+
    scale_fill_continuous(type = "viridis") + theme_bw()+ ggtitle("Core essential genes\n normal approximation")
  pNE2 <- ggplot(dfNE2, aes(x=samplesmean_NE,y=samplessd_NE) ) +theme(axis.text.x = element_text(size=13),
                                                                      axis.text.y = element_text(size=13),axis.title.x = element_text(size=12),
                                                                      axis.title.y=element_text(size=12))+
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  + xlab("mean") + ylab("sd")+
    scale_fill_continuous(type = "viridis") + theme_bw()+ ggtitle("Nonessential genes\n normal approximation")
  pFull1 <- ggarrange(pE1,pNE1,ncol=2, nrow=1, common.legend = TRUE, legend="right")
  pFull2 <- ggarrange(pE2,pNE2,ncol=2, nrow=1, common.legend = TRUE, legend="right")
  pFull <- ggarrange(pFull1,pFull2,ncol=1,nrow=2,common.legend=F)
  pdf(file=paste(plot_file_name,"_priors_data.pdf",sep=""),width=unit(6,"cm"),height=unit(5,"cm"))
  print(pFull,dpi=1200)
  dev.off()}
  return(output)
}

#' @import ggplot2
#' @importFrom rstan sampling
#' @importFrom rstan ess_tail
#' @importFrom rstan Rhat
#' @importFrom rstan ess_bulk
#' @export
core_function_1gene <- function(lfc_gene, gene_name,prior,keep_samples=F,plot_posterior=F,file_name=NA, min_tail_ESS=500)
{
  expensiveSampling=F
  veryExpensiveSampling=F
  SigmaEssential <- matrix(prior$essential$par[9:12],nrow=2)
  SigmaNE <- matrix(prior$NE$par[9:12],nrow=2)
  meanEssential <- prior$essential$par[5:6]
  meanNE <- prior$NE$par[5:6]
  samples <- sampling(stanmodels$BSure,data = list(M = length(as.vector(lfc_gene)), y=as.vector(lfc_gene),
                                           meanEssential=meanEssential,meanNE=meanNE,SigmaNE=SigmaNE,SigmaEssential=SigmaEssential),chains=2, cores=1,
                      iter=8000,thin=2,control = list(adapt_delta=0.9),refresh=8000)
  A <- as.matrix(samples)
  essTail <- apply(A,2,ess_tail)
  rhatt <- apply(A,2,Rhat)
  print(c(essTail,rhatt))
  if (max(rhatt) > 1.1 || min(essTail) < min_tail_ESS)
  {
    expensiveSampling = T
    samples <- sampling(stanmodels$BSure,data = list(M = length(lfc_gene), y=as.vector(lfc_gene),
                                             meanEssential=meanEssential,meanNE=meanNE,SigmaNE=SigmaNE,SigmaEssential=SigmaEssential)
                        ,chains=2, cores=1,iter=40000,thin=10,control = list(adapt_delta=0.999),refresh=80000)
  }
  if (max(rhatt) > 1.1 || min(essTail) < min_tail_ESS)
  {
    veryExpensiveSampling = T
    samples <- sampling(stanmodels$BSure,data = list(M = length(lfc_gene), y=as.vector(lfc_gene),
                                             meanEssential=meanEssential,meanNE=meanNE,SigmaNE=SigmaNE,SigmaEssential=SigmaEssential)
                        ,chains=2, cores=1,iter=240000,thin=60,control = list(adapt_delta=0.999),refresh=80000)
  }
  output <- c()
  posterior <- as.matrix(samples)
  output$ess_tail <- apply(posterior,2,ess_tail)
  output$ess_bulk<- apply(posterior,2,ess_bulk)
  output$rhat <- apply(posterior,2,Rhat)
  print(output$rhat)
  print(output$essTail)
  output$mean <- apply(posterior,2,mean)
  output$quant05 <- apply(posterior,2,function(x) return(quantile(x,0.05)))
  output$quant025 <- apply(posterior,2,function(x) return(quantile(x,0.025)))
  output$quant95 <- apply(posterior,2,function(x) return(quantile(x,0.95)))
  output$quant975 <- apply(posterior,2,function(x) return(quantile(x,0.975)))
  output$gene_names <- gene_name
  output$expensive_sampling <- expensiveSampling
  output$very_expensive_sampling <- veryExpensiveSampling
  if (keep_samples == T)
  {
    output$posterior <- posterior
  }
  if (plot_posterior == T)
  {
    e <- posterior[,1]
    df <- data.frame(e=e)
    pdf(file=file_name,height=8.8/2.54,width = 8.8/2.54)
    pp <- ggplot(df) + stat_density(aes(x=e),color="blue",geom="line",show.legend = T)  + xlab("essentiality score")+
      theme_classic(base_size=15)+ scale_y_continuous(expand = c(0, 0))+
      geom_vline(xintercept=as.vector(as.vector(lfc_gene)),color="grey",show.legend = T)+
      theme(legend.title=element_blank())+theme(legend.position="top")+ggtitle(gene_name)
    print(pp)
    dev.off()
  }
  return(output)
}

#' Run the BSure algorithm, which samples for each gene from the posterior distribution of gene essentiality and standard deviation
#' parameter. Furthermore it provides credible intervals of the parameters, as well as convergence checks and information
#' on whether (for each gene separately) a more or less expensive algorithm was required to achieve convergence.
#' The output is saved as an .rda file and then used by postprocessing functions.
#' @title runBSure
#' @param lfc Mx(R+2) matrix or data frame of gRNA log-fold changes, where M is the number of gRNAs and R the number of replicates,
#' the first column of the matrix contains the names of the gRNAs, the second gene names or identifiers
#' corresponding to the gRNAs in the first column
#' @param save_file_name name of the file (without .rda extension) to save the output
#' @param min_tail_ESS minimum cutoff for the tail essential sample size, to decide on whether to chose a more
#' expensive sampling algorithm, the default is set to 2000, the recommended minimum is 500. To increase the
#' speed of the algorithm reduce from the default to a minimum of 500.
#' @param vector_of_genes optional vector of genes to apply the BSure to only a subset of genes. Default is
#' NULL, in which case the algorithm is applied to all genes.
#' @param plot_folder_name If this parameter is not set to default NULL, a folder named plot_folder_name is created,
#' plots to illustrate the posterior distriubtion of the essentiality score for each gene are stored in the folder.
#' It is recommended that this option is only used with a limited number of genes, as specified by vector_of_genes.
#' @return The output is saved as an .rda file. It contains a large amount of information that is used for
#' plotting, assessment of screen quality, etc. by the postprocessing functions. For each gene it contains
#' ess_bulk bulk essential sample size, number of independent samples corresponding to the total of the
#' (dependent) samples drawn using no-U-turn-sampling
#' rhat rhat convergence criterion (Brooks and Gelman, 1998); we recommend that this should be below 1.1
#' ess_tail tail essential sample size, like bulk essential sample size, but giving information on how
#' well the algorithm was able to sample from the tails of the posterior distribution
#' mean mean of the posterior distributions of gene essentiality, standard deviation parameter and log-likelihood
#' quant05 0.05 quantiles of the posterior distributions
#' quant025 0.025 quantiles of the posteirior distributions
#' quant95 0.95 quantiles of the posterior distributions
#' quant975 0.975 quantiles of the posterior distributions
#' gene_names names of the genes in order corresponding to the output
#' expensive_sampling binary variable for each gene, indicating whether more a more expensive sampling algorithm
#' was needed
#' very_expensive_sampling binary variable for each gene, indicating whether more the most expensive sampling algorithm
#' was needed
#' probability_essential_II probability of the posterior distribution of the essentiality being
#' shifted less than 1/3 compared to a typical core-essential gene
#' probability_essential_I probability of a gene being not nonessential (probability of not being nonessential)
#' @references Brooks S, Gelman A.  General methods for monitoring convergence of iterative simulations.
#' \empph{J Comput Graph Stat.} 1998;7(4)434–455. doi10.1080/10618600.1998.10474787
#' @author Magdalena Strauss
#' @import parallel
#' @examples data(HT29_lfc_small)
#' runBSure(lfc_small,save_file_name = "temp",n_cores = 2,min_tail_ESS = 500)
#' @export
runBSure <- function(lfc,save_file_name,n_cores,min_tail_ESS=500,vector_of_genes=NULL,
                     plot_folder_name=NULL)
{
  gRNAToGene <- as.matrix(lfc[,1:2])
  lfc <- lfc[,-(1:2),drop=F]
  M <- matrix(nrow=nrow(lfc),ncol=ncol(lfc))
  M[] <- as.numeric(as.matrix(lfc))
  lfc <- M
  geneNames <- unique(gRNAToGene[,2])
  data(auxData)
  coreEssGenes <- intersect(CEGv2,coreFitnessBehanEtAl)
  nonEssGenes <- NEG
  geneNames <- unique(gRNAToGene[,2])
  x1 <- which(gRNAToGene[,2]%in%coreEssGenes)
  x2 <- which(gRNAToGene[,2]%in%nonEssGenes)
  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]
  gRNAEssentialMedian <- apply(lfc[x1,,drop=F],2,function(x) return(median(x,na.rm=T)))
  gRNA_NE_median <- apply(lfc[x2,,drop=F],2,function(x) return(median(x,na.rm=T)))
  a <- -1
  b <- 0
  k <- (a-b)/(gRNAEssentialMedian-gRNA_NE_median)
  d <- a - k*gRNAEssentialMedian
  lfc_medNorm <- matrix(0,nrow=nrow(lfc),ncol=ncol(lfc))
  lfc <- as.matrix(lfc)
  for (j in 1:ncol(lfc))
  {
    lfc_medNorm[,j] <- lfc[,j]*k[j]+d[j]
  }
  dim(lfc_medNorm) <- dim(lfc)
  lfc <- lfc_medNorm
  meanSdsEssential <- compute_mean_lfc(lfc[x1,,drop=F],gRNAToGene[x1,2])
  meanSdsNE <- compute_mean_lfc(lfc[x2,,drop=F],gRNAToGene[x2,2])
  priorDists <- fit_bivariate_nv2(meanSdsEssential$mean_lfc,meanSdsNE$mean_lfc,meanSdsEssential$sds_lfc,
                                meanSdsNE$sds_lfc,save_file_name)
  runInfStore300Samples <- function(geneName)
  {
    if (any(gRNAToGene[,2]==geneName))
    { lfc_gene <- as.vector(lfc[gRNAToGene[,2]==geneName,,drop=F])
    lfc_gene <- lfc_gene[!(is.na(lfc_gene))]
    output <- core_function_1gene(lfc=lfc_gene, geneName, prior=priorDists,keep_samples=T,
                                plot_posterior=F,file_name= paste(save_file_name,".pdf",sep=""),min_tail_ESS=min_tail_ESS)
    output$posterior <- output$posterior[seq(1,nrow(output$posterior),by = floor(nrow(output$posterior)/300)),1]
    return(output)}
  }
  mcluster <- parallel::makeCluster(n_cores)
  parallel::clusterExport(mcluster, c("runInfStore300Samples","core_function_1gene","quantile","lfc","priorDists","geneNames",
                            "sampling","ess_tail","ess_bulk","Rhat","gRNAToGene","geneNames","save_file_name",
                            "coreEssGenes","nonEssGenes","min_tail_ESS"),envir=environment())

  outputNE <- parSapply(mcluster,intersect(geneNames,nonEssGenes),runInfStore300Samples)
  save(outputNE,file= paste(save_file_name,"_EssentialNE.rda",sep=""))
  outputEss <- parSapply(mcluster,intersect(geneNames,coreEssGenes),runInfStore300Samples)
  save(outputEss,file= paste(save_file_name,"_EssentialEss.rda",sep=""))
  distNE <- c()
  for (j in 1:length(intersect(geneNames,nonEssGenes)))
  {
    distNE <- cbind(distNE,outputNE[,j]$posterior)
  }
  distEss <- c()
  for (j in 1:length(intersect(geneNames,coreEssGenes)))
  {
    distEss <- cbind(distEss,outputEss[,j]$posterior)
  }

  save(distNE,distEss,file=paste(save_file_name,"distEssNE.rda",sep=""))
  plot_posterior <- F
  if (!(is.null(plot_folder_name)))
  {
    dir.create(plot_folder_name)
    plot_posterior <- T
  }
  parallel::clusterExport(mcluster, c("distNE","distEss","plot_folder_name","plot_posterior"),envir=environment())
  runInfCompare <- function(geneName)
  {
    if (any(gRNAToGene[,2]==geneName))
    {
      lfc_gene <- as.vector(lfc[gRNAToGene[,2]==geneName,,drop=F])
      lfc_gene <- lfc_gene[!(is.na(lfc_gene))]
      save_file_name <- paste0(plot_folder_name,"/",geneName)
      output <- core_function_1gene(lfc=lfc_gene, geneName, prior=priorDists,keep_samples=T,
                                  plot_posterior=plot_posterior,file_name= paste0(save_file_name,".pdf"),min_tail_ESS=min_tail_ESS)
      lNE <- floor(0.5*min(length(distNE),length(output$posterior[,1])))
      xy <- sample(1:length(distNE),2*lNE,F)
      xx <- xy[1:lNE]
      yy <- xy[(lNE+1):(2*lNE)]
      distNE1 <- distNE[xx]
      distNE2 <- distNE[yy]
      lE <- floor(0.5*min(length(distEss),length(output$posterior[,1])))
      xy1 <- sample(1:length(distEss),2*lE,F)
      distEss <- distEss[xy1]
      samplesOutputEss <- sample(output$posterior[,1],2*lE,replace=F)
      samplesOutputNE <- sample(output$posterior[,1],lNE,replace=F)
      output$probability_essential_II <- sum(sort(samplesOutputEss) < sort(distEss) + 1/3)/length(distEss)
      output$probability_essential_I <- sum(samplesOutputEss[1:lE] < distNE1 & distNE2 > distNE1)/sum(distNE2>distNE1)
      if (!(is.null(plot_folder_name)))
      {
        pdf(paste0(save_file_name,".pdf"),height=8.8/2.54,width = 13/2.54)
        y <- c(distEss,distNE,output$posterior[,1])
        ind <- factor(c(rep("core essential",length(distEss)),rep("nonessential",length(distNE)),
                 rep(geneName,nrow(output$posterior))),levels=c(geneName,"core essential","nonessential"))
        df <- data.frame(y=y,ind=ind)
        pp <- ggplot(df,aes(x=y)) + geom_density(aes(color=ind))+theme_classic(base_size=15)+ scale_y_continuous(expand = c(0, 0))+
          theme(legend.title=element_blank())+theme(legend.position="top")+theme(plot.title = element_text(size=13))+
          ggtitle(paste0("essential I: ",toString(round(output$probability_essential_I,2)),", essential II: ",
                         toString(round(output$probability_essential_II,2))))
        print(pp)
        dev.off()
      }
      output$posterior <- NULL
      return(output)}
  }
  if (is.null(vector_of_genes))
  {vector_of_genes <- geneNames}
  output <- parSapply(mcluster,vector_of_genes,runInfCompare)
  save(output,file= paste(save_file_name,".rda",sep=""))
  stopCluster(mcluster)
}

