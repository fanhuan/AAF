# This R function estimates genome size from kmer_diversity (Fan et al 2015)
genomesize <- function(total_kmer, kmer_diversity,k,read_length,total_bp) {
  kmer_cov <- total_kmer/kmer_diversity
  base_cov <- kmer_cov * read_length / (read_length - k + 1)
  gsize <- total_bp/base_cov
  return(gsize)
}