#' Create a long file from a wide file
#'
#' Uses a minimal number  of arguments to create a long file using \code{stats::\link{reshape}}. Produces output even when long variable names and time values are not fully crossed.
#'
#' \code{tolong} is intended for the simple case in which 'wide' variables in the input data frame are identified
#' by the fact that they contain a separator character that separates the name of the variable in the long file from
#' the value of the 'time' variable that identifies the corresponding row in the long file, e.g \code{x_1, x_2, x_3} or
#' \code{brain.volume_left, brain.volume_right}.  If the separator ('_' by default) occurs in other variables, it must
#' be temporarily substituted.
#'
#' \code{\link{reshape}} does not work if long variable names and time values are not fully crossed, e.g \code{x_1, x_2, x_3, y_1, y_2}. By default \code{long}
#' creates additional variables with "NAs" so the set of variables given to \code{\link{reshape}} is fully crossed, e.g. adding a variable \code{y_3 <- NA}.
#'
#' @param data wide data frame
#' @param sep (default '_') single character separator between long names and 'time' value.
#'        Variable names with this separator are transformed to long variables.
#' @param timevar (default 'time') names of variable in output long file to identify occasions. Its values are taken from the suffix following
#'        the 'sep' character in each time-varying variable.
#' @param idvar  (default: 'id') the variable name used in the output long file to identify rows in the input wide file. It may exist in the input wide file and must, in that case,
#'        have a unique value in each row. If it does not exist, it will be created with values equal to the row numbers of the input wide file.
#' @param ids  (default \code{1:nrow(data)}) values for idvar in long file if the variable \code{idvar} does not exist in in the input wide file. Ignored
#'        if \code{idvar} exists in \code{data}.
#' @param expand (default TRUE): if 'time' values are inconsistent, fill in missing 'time's with NAs.
#' @param safe_sep temporary safe? separator
#' @param reverse (default FALSE) if TRUE, the 'time' value precedes the 
#'        variable name
#' @param ... additional parameters are passed to \code{\link{reshape}}.
#' @return 'long' file with each wide row repeated as many times as there are distinct values for the 'timevar' variable.
#' @seealso \code{\link{towide}} for many examples using both 'towide' and 'tolong'.
#' @examples
#' z <- data.frame(id =letters[1:10], id2= 11:20, v_L = 1:10, v_R = 11:20)
#' tolong(z)
#' tolong(z, timevar = 'Side', idvar = 'idn', ids = LETTERS[1:10])
#' tolong(z, timevar = 'Side', idvar = 'idn', ids = z$id2)
#'
#' # unbalanced times
#' z <- data.frame(id =letters[1:10], id2= 11:20, v_L = 1:10, v_R = 11:20, z_L = 21:30)
#' tolong(z)
#'
#' # multi-character sep
#' z <- data.frame(id =letters[1:10], id2= 11:20, HPC_head_R = 1:10, HPC_tail_R = 11:20, HPC_head_L = 21:30, HPC_tail_L = 31:40)
#' names(z) <- sub("(_[LR]$)","_\\1", names(z))
#' names(z)
#' (zz <- tolong(z, sep = "__", timevar = "Side"))
#' zz$id3 <- rownames(zz)
#' tolong(zz, idvar = 'id3' ,timevar = 'Part')
#'
#' dd <- data.frame( y.a = 1:3, y.b = 1:3, x.a= 1:3, time = 1:3,
#'     x.b = 11:13, x.c = 21:23, id = c('a','a','b'))
#' tolong(dd, sep = '.')
#' dl <- tolong(dd, sep = '.', timevar = "type", idvar = 'patient')
#' towide(dl, idvar = 'patient', timevar = 'type')
#'
#' # Long file with additional constants
#'
#' dl <- data.frame(name = rep(c('A','B','C'), c(3,3,2)),
#'                  site = c('head','neck','jaw','chest')[
#'                    c(1,2,3,1,2,3,1,4)],
#'                  sex = rep(c('male','female','male'), c(3,3,2)),
#'                  var1 = 1:8,
#'                  var2 = 11:18,
#'                  invar = rep(1:3, c(3,3,2)))
#' towide(dl, c('name'), 'site')
#' #
#' # Two indexing variable: e.g. hippocampal volume 2 sides x 3 sites
#' #
#' dl <- data.frame(name = rep(LETTERS[1:3], each = 6),
#'                  side = rep(c('left','right'), 9),
#'                  site = rep(rep(c('head','body','tail'),each = 2),3),
#'                  volume = 1:18,
#'                  sex = rep(c('female','male','female'), each = 6),
#'                  age = rep(c(25, 43, 69), each = 6))
#' dl
#' (dlsite <- towide(dl, c('name','side'), 'site'))
#' (dlsite.side <- towide(dlsite, c('name'), 'side'))
#' #
#' # Flipping a data frame
#' #
#' z <- data.frame(vname = 
#'    rep(c('v1','v2','v3'), each = 4),
#'    country = rep(c('Angola','Benin','Chad','Denmark'), 3),
#'    code = rep(c('ANG','BEN','CHA','DEN'),3),
#'    val__2011 = 2011 + seq(.01,.12,.01),
#'    val__2012 = 2012 + seq(.01,.12,.01),
#'    val__2013 = 2013 + seq(.01,.12,.01),
#'    val__2014 = 2014 + seq(.01,.12,.01),
#'    val__2015 = 2015 + seq(.01,.12,.01)
#' )
#' z
#' z %>% 
#'   tolong(sep= '__') 
#' 
#' z %>% 
#'   tolong(sep= '__', timevar = 'year') %>% 
#'   .[!grepl('^id$',names(.))] %>% 
#'   towide(timevar = 'vname', idvar = c('code','year'))
#' 
#' \dontrun{
#' # Extracting chains from a stanfit object in the 'rstan' package
#' # If 'mod' is a stanfit model
#' library(rstan)
#' library(spida2)
#' df <- as.data.frame(extract(mod, permute = F))
#' dl <- tolong(df, sep = ':', reverse = T)
#' }
#' @export
tolong <- function (data, sep = "_",  timevar = 'time',
                  idvar = 'id', ids = 1:nrow(data),
                  expand = TRUE, safe_sep = "#%@!", 
                  reverse = F, ...) {
  if (timevar %in% names(data)) warning(paste("Variable",timevar, "in data is replaced by a variable to mark occasions"))
  if (idvar %in% names(data)) {
    idwide <- data[[idvar]]
    if( length(unique(idwide)) != length(idwide)) {
      warning(paste ("idvar:", idvar, "must have unique values. It will be replaced in the output data frame with a variable containing row numbers of the input data frame"))
      data[[idvar]] <- ids
    }
  }
  Flip <- function(nams, sep) {
    expr <- paste0('^(.*)(\\Q',sep,'\\E)(.*)$')
    new <- "\\3\\2\\1"
    sub(expr,new,nams)
  }
  if (reverse) names(data) <- Flip(names(data), sep)
  if (expand) {
    namessafe <- sub(sep, safe_sep, names(data), fixed = TRUE)
    varnames <- grep(safe_sep, namessafe, value = TRUE, fixed = TRUE)
    names <- unique(sub(paste(safe_sep, ".*$", sep = ""),
                        "", varnames))
    times <- unique(sub(paste("^.*", safe_sep, sep = ""),
                        "", varnames))
    allnames <- paste(rep(names, each = length(times)), sep,
                      rep(times, length(names)), sep = "")
    z <- data
    for (nn in allnames) {
      z[[nn]] <- if (is.null(data[[nn]]))
        NA
      else data[[nn]]
    }
    data <- z
  }
  namessafe <- sub(sep, safe_sep, names(data), fixed = TRUE)
  namestimes <- sub(paste("^.*", safe_sep, sep = ""), "ZZZZZ",
                    namessafe)
  ord <- order(namestimes)
  data <- data[, ord]
  stats::reshape(data, direction = "long", sep = sep,
          varying = grep(sep, names(data), fixed = TRUE),
          timevar = timevar, idvar = idvar,
          ids = ids, ...)
}
#' @export
#' @rdname tolong
long <- tolong
#  z <- data.frame(id =letters[1:10], id2= 11:20, v_L = 1:10, v_R = 11:20)
#  long(z)
#  long(z, timevar = 'Side', idvar = 'idn', ids = LETTERS[1:10])
#  long(z, timevar = 'Side', idvar = 'idn', ids = z$id2)
# z$nid <- 1:5
# long(z, idvar = 'nid')
# z
#  # unbalanced times
#  z <- data.frame(id =letters[1:10], id2= 11:20, v_L = 1:10, v_R = 11:20, z_L = 21:30)
#  long(z,idvar='rowid')
#
#  reshape(direction= 'long', z, sep = '_', v.names= 'v') # does not work
