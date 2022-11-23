#' gut microbe OTU data (species level)
#'
#'The original nucleotide sequences of this study were deposited to the NCBI
#'Sequence Read Archive under accession number SRP128619.
#'@format A data frame with 65 rows and 21 column, contain first column as microbe ID:
#'@docType data
#'@name gut_microbe
#'@usage data(gut_microbe)
"gut_microbe"

#' mustard microbe OTU data
#'
#'Wagner, M. R. et al. Host genotype and age shape the leaf and root microbiomes of a wild perennial plant.
#'Nat. Commun. 7:12151 doi: 10.1038/ncomms12151 (2016)
#'This dataset is a subset of otuTable97, we select location = JAM, keep samples with both root
#'and leaf data, and then run data_cleaning first (set x = 50) to reduce size of this data.
#'Moreover, sample 8_1382 is removed for the outlier reason.
#'@format A data frame with 1557 rows and 176 column, contain first column as OTU ID:
#'@docType data
#'@name mustard_microbe
#'@usage data(mustard_microbe)
"mustard_microbe"
