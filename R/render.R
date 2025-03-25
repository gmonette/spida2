#' Render with Rmarkdown and keep intermediate files
#' 
#' Following https://github.com/rstudio/rmarkdown/issues/17 
#' add the YAML: \code{knit: spida2::render_keep} to a .R Rmarkdown 
#' document with output YAML \code{output: html_document}.
#' 
#' @param input file to be rendered
#' @param ... other arguments passed to \code{\link{rmarkdown::render}}.
#' @export
render_keep <- function(input, ...) {
  rmarkdown::render(input, clean = FALSE, envir = new.env(), ...)
}