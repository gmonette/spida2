#' Create a long file from a wide file
#'
#' Uses a minimal number of arguments to create a long file
#' using \code{stats::\link{reshape}}. Produces output even
#' when long variable names and time values are not fully
#' crossed.
#'
#' \code{tolong} is intended for the simple case in which
#' 'wide' variables in the input data frame are separator string that
#' separates the name of the variable in the long file from
#' the value of the 'time' variable that identifies the
#' corresponding row in the long file, e.g \code{x_1, x_2,
#' x_3} or \code{brain.volume_left, brain.volume_right}.  Since
#' the separator ('_' by default) may occur in other variables,
#' \code{tolong} offers two mechanisms to avoid misinterpreting
#' those occurrences as separators. If there are multiple
#' occurrences of the separator string in a variable name, only the 
#' last occurrence is interpreted as a separator.  Secondly,
#' the \code{valuepattern} parameter can specify a regular expression
#' to identify allowed 'time' value strings at the end of variable 
#' names.  The common case where the 'time' values are numerical
#' can be specified with the \code{numericalpattern} parameter.
#'
#' \code{\link{reshape}} does not work if long variable names
#' and time values are not fully crossed, e.g \code{x_1, x_2,
#' x_3, y_1, y_2}. By default \code{long} creates additional
#' variables with "NAs" so the set of variables given to
#' \code{\link{reshape}} is fully crossed, e.g. adding a
#' variable \code{y_3 <- NA}.
#'
#' Compare the functionality of \code{tolong} with that of
#' \code{\link{tidyr::gather}} and of
#' \code{\link{tidyr::pivot_longer}}. 'tolong' depends on the
#' format of variable names to identify variables whose values
#' become new variables in the long form of the data and which
#' labels are used as the indices of the indexing variable,
#' whose default name is 'time', which can be set to another
#' value with the "timevar" argument. "tolong" can handle many
#' 'time-varying' variables. "gather" can only handle one.
#' "pivot_longer" can handle many and might be considered a
#' replacement for "to_long" which has the disadvantage of
#' frequently requiring the renaming of variables, an easier
#' task for those who have mastered the use of regular
#' expressions, but potentially challenging otherwise.
#' 
#' @param data wide data frame
#' @param sep (default '_') single or multiple character separator between
#'   long names and 'time' value. Variable names with this
#'   separator are transformed to long variables. If the string
#'   occurs multiple times in a variable name, only the last
#'   occurrence is treated as a separator.
#' @param valuepattern a regular expression to match the form of time values
#'   at the end of variable names immediately following
#'   the separator.  Specifying this pattern can avoid
#'   misinterpreting separators in variable names that
#'   are not intended to turned into long variables.
#' @param numericalpattern (default FALSE) if TRUE, 
#'   \code{valuepattern} is set to "[0-9]+$" to match
#'   numerical values as time values at the end of
#'   variable names to be turned into long variables. 
#'   This parameter is only
#'   necessary if the separator occurs in variable names where
#'   it is not intended to identify a long variable. 
#' @param timevar (default 'time') names of variable in output
#'   long file to identify occasions. Its values are taken
#'   from the suffix following the 'sep' character in each
#'   time-varying variable.
#' @param idvar  (default: 'id') the variable name used in the
#'   output long file to identify the provenance row in the 
#'   input wide file.
#'   It may exist in the input wide file and must, in that
#'   case, have a unique value in each row. If it does not
#'   exist, it will be created with values equal to the row
#'   numbers of the input wide file.
#' @param ids  (default \code{1:nrow(data)}) values for idvar
#'   in long file if the variable \code{idvar} does not exist
#'   in the input wide file. Ignored if \code{idvar} exists
#'   in \code{data}.
#' @param expand (default TRUE): if 'time' values are
#'   inconsistent, fill in missing 'time's with NAs.
#' @param safe_sep temporary safe? separator
#' @param reverse (default FALSE) if TRUE, the 'time' value
#'   precedes the variable name
#' @param ... additional parameters are passed to
#'   \code{\link{reshape}}.
#' @return 'long' file with each wide row repeated as many
#'   times as there are distinct values for the 'timevar'
#'   variable. The rownames show the provenance of each row
#'   by combining the value of `id` with the value of `time`
#'   separated by a period
#' @seealso \code{\link{towide}} for many examples using both
#'   'towide' and 'tolong'.
#' @examples
#' z <- data.frame(id =letters[1:10], id2= 11:20, v_L = 1:10, v_R = 11:20)
#' z
#' tolong(z)
#' tolong(z, timevar = 'Side', idvar = 'idn', ids = LETTERS[1:10])
#' tolong(z, timevar = 'Side', idvar = 'idn', ids = z$id2)
#'
#' # unbalanced times
#' z <- data.frame(id =letters[1:10], id2= 11:20, v_L = 1:10, v_R = 11:20, z_L = 21:30)
#' z
#' tolong(z)
#'
#' # a separator with multiple occurrences:
#' z <- data.frame(id =letters[1:10], id2= 11:20, v_a_L = 1:10, v_a_R = 11:20, z_L = 21:30)
#' z
#' # The previous version of tolong() would have produced an error 
#' # due to multiple occurrences of the default separator '_'
#' # but the new version matches only the last occurrence in
#' # each variable name. The sublast() function helps by
#' # matching only the last occurrence to facilitate 
#' # replacing it with a new unique separator, but it is
#' # no longer necessary to do this:
#' zz <- z
#' names(zz) <- sublast('_', '__', names(zz))
#' tolong(zz, sep = '__')
#' # or, now,, with the same result:
#' tolong(z)
#' #
#' # - sep can use more than one character
#' # - the character string is interpreted literally, 
#' # i.e. if special regular expression characters
#' # they are interpreted literally.
#' z <- data.frame(id =letters[1:10], id2= 11:20, HPC_head_R = 1:10, HPC_tail_R = 11:20, HPC_head_L = 21:30, HPC_tail_L = 31:40)
#' z
#' names(z) <- sub("(_[LR]$)","_\\1", names(z))
#' names(z)
#' (zz <- tolong(z, sep = "__", timevar = "Side"))
#' zz$id3 <- rownames(zz)
#' tolong(zz, idvar = 'id3' ,timevar = 'Part')
#'
#' dd <- data.frame( y.a = 1:3, y.b = 1:3, x.a= 1:3, time = 1:3,
#'     x.b = 11:13, x.c = 21:23, id = c('a','a','b'))
#' dd
#' tolong(dd, sep = '.')
#' dl <- tolong(dd, sep = '.', timevar = "type", idvar = 'patient')
#' dl
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
#' dd <- data.frame(var__tag_1 = 1:10, var__tag_2 = 1:10)
#' dd <- data.frame(var__tag_1 = 1:10, var__tag_2 = 1:10)
#' dd <- data.frame(var.1 = 1:10, var.2 = 1:10)
#' tolong(dd, sep = '__')
#' tolong(dd, sep = '.')
#' 
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
                    valuepattern = if(numericalpattern) '[0-9]+' else '.*',
                    numericalpattern = FALSE,
                    expand = TRUE, safe_sep = '@3-2861-2579@', 
                    reverse = F, ...) {
  data <- as.data.frame(data) # if data is a tibble then reshape fails
  if (timevar %in% names(data)) warning(paste("Variable",timevar, "in data is replaced by a variable to mark occasions. Use the 'timevar' argument to specify a different variable name in the output data frame"))
  if (idvar %in% names(data)) {
    idwide <- data[[idvar]]
    if( length(unique(idwide)) != length(idwide)) {
      warning(paste ("idvar:", idvar, "must have unique values. It will be replaced in the output data frame with a variable containing row numbers of the input data frame. Use the 'idvar' argument to specify a different variable name in the output data frame"))
      data[[idvar]] <- ids
    }
  }
  
  # sample(100000000,3)
  
  # s1 <- '@1-7418-1770@'
  # s2 <- '@2-1549-5289@'
  # s3 <- safe_sep
  # 
  # n0 <- names(data)
  # n1 <- gsub(sep, s1, n0, fixed = TRUE)  # sep to s1  (need to protect if _)
  # n2 <- gsub('_', s2, n1)                # '_' if not sep to s2
  # n3 <- gsub(s1, '_', n2)                # sep to '_' which has been protected
  # n4 <- sub(paste0('(_)(?!.*_.*)(',valuepattern,')$'),paste0(s3,"\\2"), n3, perl = TRUE)  # last sep to s3 + final string
  # n5 <- gsub('_',safe_sep, n4)                # remaining '_'s back to sep     
  # 
  
  str1 <- '@-3480-3285-@'
  # str2 <- '@-0752-8427-@'
  str2 <- safe_sep
  
  # sep <- '__'
  # valuepattern <- '.*'
  # (n0 <- c('val__a.1','val__a.2','val__a.3','val__a.4'))
  # (n0 <- c('val__a.1','val__a.2','val__a.a','val__a.b'))
  n0 <- names(data)
  n1 <- gsub(sep, str1, n0, fixed = TRUE)
  reppat <- paste0('(',str1,')(?!.*',str1,'.*)(',valuepattern,')$')
  n2 <- sub(reppat, paste0(str2,"\\2"), n1, perl = TRUE)  # last sep to s3 + final string
  n3 <- gsub(str1, sep, n2)
  
  names(data) <- n3
  if(expand){
    # create variables with NAs for missing time values
    # namessafe <- sub(sep, safe_sep, names(data), fixed = TRUE)
    varnames <- grep(safe_sep, names(data), value = TRUE, fixed = TRUE)  # pick vars with safe_sep
    names <- unique(sub(paste(safe_sep, ".*$", sep = ""),                # get root names
                        "", varnames, fixed = FALSE))
    
    times <- unique(sub(paste("^.*", safe_sep, sep = ""),                # get values
                        "", varnames, fixed = FALSE))
    allnames <- paste(rep(names, each = length(times)), safe_sep,        # generate all names
                      rep(times, length(names)), sep = "")
    newnames <- setdiff(allnames, varnames)
    
    for (nn in newnames) data[[nn]] <- NA
  }
  
  ret <- stats::reshape(data, direction = "long", sep = safe_sep,
                        varying = grep(safe_sep, names(data), fixed = TRUE),
                        timevar = timevar, idvar = idvar,
                        ids = ids,...)
  # names(ret) <- gsub(s2, '_', names(ret), fixed = TRUE)
  ret[order(ret[[idvar]]),]
}

if(FALSE){
  # test tolong_new
  # 
  # Test finding last separator when separator not _
  (zd <- data.frame(a=1:4, b.left.1 = 1:4, b.left.2 =1:4, b.right.1 = 1:4, b.right =1:4 , time = 1:4))
  tolong(zd, sep = '.')  # wrong
  tolong(zd, sep = '.', val = "[0-9]+")  # correct
  tolong(tolong(zd, sep = '.', val = "[0-9]+"),sep = '.', timevar = 'side')  # correct
  
  
  zd <- data.frame(a=1:4, b__1 = 1:4, b__2=1:4, c__1 = 1:4, time = 1:4)
  tolong(zd, sep = '__')
  tolong(zd, sep = '__', timevar = 'index')
  tolong(zd, sep = '__', timevar = 'occasion', idvar = 'row')
  zd <- data.frame(a=1:4, b_left_1 = 1:4, b_right_1=1:4,b_left_2 = 1:4, b_right_2=1:4, c_1 = 1:4, time = 1:4)
  (zdd <- tolong(zd, sep = '_'))
  (tolong(zdd, sep = '_'))
  (tolong(zdd, sep = '_', timevar = 'side', idvar = 'row'))
  zd <- data.frame(a=1:4, b_left_1 = 1:4, b_right_1=1:4,b_left_2 = 1:4,  c_1 = 1:4, time = 1:4)
  (zdd <- tolong(zd, sep = '_'))
  (tolong_new(zdd, sep = '_'))
  zd <- data.frame(a=1:4, b.left.1 = 1:4, b.right.1=1:4,b.left.2 = 1:4,  c.1 = 1:4, time = 1:4)
  (zdd <- tolong(zd, sep = '.'))
  (tolong(zdd, sep = '.'))
  
}
#' @export
#' @rdname tolong
tolong_old <- function (data, sep = "_",  timevar = 'time',
                        idvar = 'id', ids = 1:nrow(data),
                        expand = TRUE, safe_sep = "#%@!", 
                        reverse = F, ...) {
  if (timevar %in% names(data)) warning(paste("Variable",timevar, "in data is replaced by a variable to mark occasions. Use the 'timevar' argument to specify a different variable name"))
  if (idvar %in% names(data)) {
    idwide <- data[[idvar]]
    if( length(unique(idwide)) != length(idwide)) {
      warning(paste ("idvar:", idvar, "must have unique values. It will be replaced in the output data frame with a variable containing row numbers of the input data frame. Use the 'idvar' argument to specify a different variable name"))
      data[[idvar]] <- ids
    }
  }
  # Check for multiple occurrences of the separator in variable names
  
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

