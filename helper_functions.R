#' Clean and prepare the data from IMNGS
#'
#' Divides the taxonomy into a new matrix for each otu
#'
#' @param taxonomy Last column of files a string ; separated with domain,
#' phylum, vlass, order, family, genus and species.
#' @param otus The name of the rows
#'
#' @return
#' A matrix with the taxonomic information ready for the package phylo
taxonomy <- function(taxonomy, otus){
    taxonomy <- sapply(as.character(taxonomy), strsplit, split = ";")
    names(taxonomy) <- otus
    otus_tax <- t(sapply(taxonomy, '[', seq(max(sapply(taxonomy, length)))))
    colnames(otus_tax) <- c("Domain", "Phylum", "Class", "Order",
                            "Family", "Genus", "Species")
    # Remove spaces
    otus_tax <- apply(otus_tax, 1:2, sub, pattern = "\\s", replacement = "")
    otus_tax <- apply(otus_tax, 1:2, sub, pattern = "[;:]", replacement = "")
    otus_tax <- apply(otus_tax, 1:2, sub, pattern = "^([a-z]__)", replacement = "")
    otus_tax[otus_tax == ""] <- NA # Remove empty cells
    otus_tax
}
