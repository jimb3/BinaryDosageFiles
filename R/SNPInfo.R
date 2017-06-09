#' Function to return SNP information read by ExtractDosages'
#' 
#' Function to return information on the SNPs read in by ExtractDosages
#' or ExtractMoreDosages
#' 
#' @param d
#' List returned from ExtractDosages or ExtractMoreDosages
#' @param x
#' Integer or vector of integers indicates SNPs to retrieve data on
#' @return 
#' Data frame with chromosome, SNPName, location in base pairs, reference allele, and alternate allele
#' @export
SNP_Info <- function(d, x) {
  Chromosome = d$Inputs$MapData$Chromosome[d$SNPID[x]]
  SNPName = d$Inputs$MapData$SNP[d$SNPID[x]]
  Location = d$Inputs$MapData$BasePairs[d$SNPID[x]]
  ReferenceAllele = d$Inputs$MapData$Allele1[d$SNPID[x]]
  AlternateAllele = d$Inputs$MapData$Allele2[d$SNPID[x]]
  return (data.frame(Chromosome,
               SNPName,
               Location,
               ReferenceAllele,
               AlternateAllele))
}