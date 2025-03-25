#' Read and print bibliographical entries
#' 
#' The following doesn't work
#' 
#' Consider packages: knitr, bibtex, vitae, etc.
#' 
#' Problem is that bibtex::read.bib rejects entries with e.g. 'online' type 
#' and entries that are incomplete which makes identifying them
#' and correcting them much more awkward.
#' 
# 
# read.bib.all <- function(bibfile, ok=c("article", "book", "booklet", "inbook", "incollection", 
#                                        "inproceedings", "manual", "mastersthesis", "misc", "phdthesis", 
#                                        "proceedings", "techreport", "unpublished")) {
#   bib <- bibtex::do_read_bib(bibfile) 
#   bib <- lapply(bib, function(x) {
#     if(!(attr(x,'entry') %in% ok)) attr(x,'entry') <- 'misc'
#     x[['bibtype']] <- attr(x,'entry')
#     x[['key']] <- attr(x,'key')
#     x
#   })
#   names(bib) <- sapply(bib, attr, 'key')
#   
#   bib
# }