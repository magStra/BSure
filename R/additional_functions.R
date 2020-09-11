#' @title create_small_subset
#' This function takes as input a matrix of count data from a CRISPR screen
#' and returns a smaller subset containing 20 core essential and 20 nonessential genes,
#' as well as 20 genes known to be differentially essential across different (cancer) cell lines
#' (Behan et al., 2019). The purpose of this function is to create a small sample dataset to illustrate
#' the method.
#' @references Behan FM et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens.
#' Nature. 2019;568(7753):511–516.
#' @param counts Mx(R+2) matrix or data frame of gRNA counts, where M is the number of gRNAs and R the number of replicates,
#' the first column of the matrix contains the names of the gRNAs, the second gene names or identifiers
#' corresponding to the gRNAs in the first column
#' @return reduced matrix of gRNA counts, containing 20 core-essential genes, 20 nonessential genes, and 20
#' genes know to be differentially essential across different (cancer) cell lines
#' @export
#' @examples data("counts_HT29")
#' counts_HT29_small_subset <- create_small_subset(counts_HT29)
create_small_subset <- function(counts)
{
  data("auxData")
  coreEssGenes <- intersect(CEGv2,coreFitnessBehanEtAl)
  nonEssGenes <- NEG
  data("Behan_et_al_list_differentially_essential_sub")
  genes <- counts[,2]
  gene_subset <- c(sample(coreEssGenes,20,replace=F),
   sample(nonEssGenes,20,replace=F),
   sample(Behan_et_al_list_differentially_essential_sub,20,replace=F))
  subset_counts <- counts[counts[,2]%in%gene_subset,]
  return(subset_counts)
}
