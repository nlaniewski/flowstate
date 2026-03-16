FCSoffsets <- function(fcs.file.path, return.version = FALSE, update.data.offsets = TRUE){
  if(grepl("fcs", tools::file_ext(fcs.file.path), ignore.case = T) & file.exists(fcs.file.path)){
    ## open a connection to the file (.fcs); read binary mode
    con <- file(fcs.file.path, open = "rb")
    on.exit(close(con))
    ## first six bytes contain "FCS#.#"; version identifier
    version <- readChar(con, 6)
    if(return.version) return(version)
    ## version test
    if(!version %in% sprintf("FCS3.%d", 0:2)){#sprintf("FCS3.%d",0:2)
      stop("Is this a .fcs file (version 3.0, 3.1, or 3.2)? flowstate has only been tested using FCS 3.0/3.1/3.2 files.")
    }
  }else{
    stop("Need a valid '.fcs' file path/file.")
  }
  ## HEADER segment
  ## version identifier; TEXT segment; DATA segment; ANALYSIS segment
  ## 6+4+8*2+8*2+8*2 = 58 bytes
  ## initialize a list for storing byte offsets
  offsets <- stats::setNames(
    vector(mode = "list", length = 4L),
    nm = c("version", "TEXT", "DATA", "ANALYSIS")
  )
  version <- substring(version, 4, nchar(version))
  offsets$version <- as.double(version)
  ## next four bytes contain space characters (ASCII 32);
  ## stop otherwise (not a legal/expected .FCS file)
  stopifnot(readChar(con, 4) == "    ")
  ## byte offsets for TEXT, DATA, and ANALYSIS
  ## character strings; 8 bytes; right justified with space padding
  ## start and end positions
  for(i in c("TEXT","DATA","ANALYSIS")){
    offsets[[i]] <- c(
      start = as.integer(readChar(con, 8)),
      end = as.integer(readChar(con, 8))
    )
  }
  ## byte offsets for OTHER segment -- if it exists
  if(offsets$TEXT[1] > 58){
    n <- offsets$TEXT[1] - 58
    seek(con,58)
    string <- readChar(con, nchars = n)
    m <- gregexpr('[0-9]+', string)
    offsets$OTHER <- as.numeric(unlist(regmatches(string, m)))
  }
  ## for DATA segments larger than 99,999,999 bytes: header offsets will be '0'
  ## update DATA offsets using keyword-value pairs stored in TEXT segment
  ## need '$BEGINDATA' and '$ENDDATA' to update offsets
  ## as attributes: '$BYTEORD' ('endianness'), '$PAR' (number of parameters)
  if(update.data.offsets){
    kv <- readFCStext(con = con, offsets = offsets[['TEXT']])
    offsets[['DATA']] <- replace(
      offsets[['DATA']],
      values = as.numeric(unlist(kv[c('$BEGINDATA', '$ENDDATA')]))
    )
    endian.string <- c("4,3,2,1" = 'big', "1,2,3,4" = 'little')
    if(!kv[['$BYTEORD']] %in% names(endian.string)){
      stop("Keyword '$BYTEORD' has unexpected endianness...data segment cannot be appropriately parsed.")
    }else{
      attr(offsets[['DATA']], "endianness") <- endian.string[[kv[['$BYTEORD']]]]
      attr(offsets[['DATA']], "par.n") <- as.numeric(kv[['$PAR']])
    }
  }
  ##
  return(offsets)
}

readFCStext <- function(fcs.file.path, return.string = FALSE, con = NULL, offsets = NULL){
  if(is.null(con) & is.null(offsets)){
    ## use the return of 'FCSoffsets(...)[["TEXT"]]'
    offsets <- FCSoffsets(fcs.file.path)[["TEXT"]]
    ## open a connection to the file (.fcs); read binary mode
    con <- file(fcs.file.path, open = "rb")
    on.exit(close(con))
  }
  ## seek the byte position where the TEXT segment starts
  seek(con, where = offsets['start'])
  ## hex/ASCII encoded keyword-value pairs
  txt <- readBin(
    con = con,
    what = "raw",
    n = diff(offsets) + 1
  )
  ## convert
  txt <- rawToChar(txt)
  ##
  if(return.string){
    return(txt)
  }
  ## delimiter
  delimiter <- substr(txt, 1, 1)
  ## split the string; start at position 2
  kv <- strsplit(substr(txt, 2, nchar(txt)), delimiter, fixed = TRUE)[[1]]#match exactly the delimiter
  ## keyword value pairs; keyword name (odd) and value (even)
  kv <- stats::setNames(
    kv[seq_along(kv) %% 2 != 1],
    nm = kv[seq_along(kv) %% 2 == 1]
  )
  ## as a list
  kv <- as.list(kv)
  return(kv)
}

add.identifier.proj <- function(keywords){
  if(any(grepl("^\\$PROJ$", names(keywords)))){
    proj <- grep("^\\$PROJ$", names(keywords), value = T)
  }else if(any(grepl("$DATE", names(keywords), fixed = T))){
    proj <- grep("$DATE", names(keywords), value = T, fixed = T)
  }
  ## return the identifier
  keywords[[proj]]
}

parameters.to.data.table <- function(
    keywords,
    add.PROJ.identifier = TRUE
)
{
  ## using the return of 'readFCStext(...)'
  ## number of parameters
  par.n <- as.numeric(keywords[['$PAR']])
  ## regex pattern to find '$P#' and 'P#'
  par.pattern <- "^\\$P\\d+|^P\\d+"
  parameter.names <- grep(
    pattern = par.pattern,
    x = names(keywords),
    value = TRUE
  )
  ## parameters; named vector
  pars <- trimws(unlist(keywords[parameter.names]))
  ## split parameters vector into 'types'; list
  par.types <- split(x = pars, f = factor(sub(par.pattern, "", names(pars))))
  ## create a data.table
  dt.parameters <- data.table::data.table(par = paste0('$P', seq(par.n)))
  ## fill the data.table by 'par.type'
  for(j in names(par.types)){
    par.vec <- par.types[[j]]
    par.i <- as.integer(gsub("\\D+", "", names(par.vec)))
    data.table::set(
      dt.parameters,
      i = par.i,
      j = j,
      value = par.vec
    )
  }
  ## add 'TYPE' keyword-value pair; encoded in CyteK/SpectroFlo files;
  ## add for other platforms; will need to be updated -- instrument-specific
  if(!'TYPE' %in% names(dt.parameters)){
    dt.parameters[grep("Time", N, ignore.case = T), TYPE := "Time"]
    dt.parameters[grep("FSC", N), TYPE := "Forward_Scatter"]
    dt.parameters[grep("SSC", N), TYPE := "Side_Scatter"]
    ## if SONY ID7000; '$CYT/LE-ID7000C';
    ## grep for '[0-9]{3}CH[0-9]+-A'
    dt.parameters[grep('[0-9]{3}CH[0-9]+-A', N), TYPE := "Raw_Fluorescence"]
    ## if CYTOF|DVS|FLUIDIGM;
    ## grep for 'Di$' ; "Center|Offset|Width|Residual" ; "Event_length"
    dt.parameters[grep('Di$', N), TYPE := "Ion_Count"]
    dt.parameters[i = grep("Center|Offset|Width|Residual", N), TYPE := "Gaussian"]
    dt.parameters[i = grep("Event_length", N), TYPE := "Event_Length"]
  }
  ## add '$PROJ' identifier; if not found, use '$DATE' instead
  if(add.PROJ.identifier){
    dt.parameters[, PROJ := as.factor(add.identifier.proj(keywords))]
  }
  ## pre-form an alias data.table for setting syntactically valid names
  ## alias allows use of NSE without having to use back ticks on unquoted variables
  alias <- data.table::data.table(
    N = dt.parameters[['N']],
    S = dt.parameters[['S']]
  )
  ## make N syntactically valid
  alias[i = !grepl("[FS]SC|Time", N), N.alias := gsub(" |-|/", "", sub("-A$", "", N))]
  alias[i = grepl("[FS]SC|Time", N), N.alias := sub("-", "_", sub("SSC-B", "SSCB", N))]
  ## make S syntactically valid
  if("S" %in% names(dt.parameters)){
    alias[, S.alias := gsub(" |-|/", "", S)]
    alias[is.na(S.alias), S.alias := N.alias]
    ##
    alias[!is.na(S), S_N.alias := paste(S.alias, N.alias, sep = "_")]
    alias[is.na(S_N.alias), S_N.alias := N.alias]
  }
  ## add as attribute to 'dt.parameters'
  data.table::setattr(dt.parameters, name = "alias", alias)
  ## return the data.table
  invisible(dt.parameters)
}

keywords.to.data.table <- function(keywords, drop.primary = TRUE, drop.spill = TRUE){
  ## using the return of 'readFCStext(...)'
  ## regex pattern to find '$P#' and 'P#' -- inverted
  par.pattern <- "^\\$P\\d+|^P\\d+"
  keyword.names<-grep(
    pattern = par.pattern,
    x = names(keywords),
    value = TRUE,
    invert = T
  )
  ## keywords; named vector
  kw <- trimws(unlist(keywords[keyword.names]))
  ## create a data.table
  dt.keywords <- data.table::as.data.table(as.list(kw))
  ## drop required keywords
  if(drop.primary){
    dt.keywords[, (fcs.text.primary.required.keywords) := NULL]
  }
  ## drop spillover keyword-value pair (string)
  if(drop.spill){
    spill.name <- grep("spill", names(dt.keywords), ignore.case = T, value = T)
    dt.keywords[, (spill.name) := NULL]
  }
  ## add keywords related to originality/modification; as defined for FCS 3.1
  dt.keywords[, '$ORIGINALITY' := "DataModified"]
  dt.keywords[, '$LAST_MODIFIED' := toupper(format(Sys.time(), "%d-%b-%Y %H:%M:%OS2"))]
  dt.keywords[, '$LAST_MODIFIER' := sprintf("flowstate_%s", utils::packageVersion("flowstate"))]
  ## return the data.table
  invisible(dt.keywords)
}

spill.to.data.table <- function(keywords, add.PROJ.identifier = TRUE){
  ## using the return of 'readFCStext(...)'
  ## spill/spillover keyword name
  spill.name <- grep('spill', names(keywords), ignore.case = TRUE, value = T)
  ##parse the spillover string
  if(length(spill.name) == 1){
    spillover.string <- keywords[[spill.name]]
    spill.split <- unlist(strsplit(spillover.string, ","))
    n.cols <- as.numeric(spill.split[1])
    col.names <- spill.split[2:(n.cols + 1)]
    vals <- as.numeric(spill.split[(n.cols + 2):length(spill.split)])
    dt.spill <- data.table::as.data.table(
      matrix(
        data = vals,
        ncol = n.cols,
        byrow = TRUE,
        dimnames = list(NULL, col.names)
      )
    )
    # rownames(dt.spill)<-names(dt.spill)
    ##add '$PROJ' identifier; if not found, use '$DATE' instead
    if(add.PROJ.identifier){
      data.table::setattr(dt.spill, 'PROJ', add.identifier.proj(keywords))
    }
    ##
    data.table::setattr(dt.spill, name = "applied", value = FALSE)
    ##
    return(dt.spill)
  }else{
    data.table::data.table()
  }
}

readFCSdata <- function(fcs.file.path){
  ## use the return of 'FCSoffsets(...)[["DATA"]]'
  offsets <- FCSoffsets(fcs.file.path)[["DATA"]]
  ## use the return of 'readFCStext(...)'
  keywords <- readFCStext(fcs.file.path)
  ## parameter names '$P#N' -- required keyword-value pair; unique name
  pars.N <- unlist(keywords[grep("\\$P\\d+N", names(keywords), value = T)])
  pars.N <- pars.N[order(as.numeric(gsub("\\D", "", names(pars.N))))]
  ## open a connection to the file (.fcs); read binary mode
  con <- file(fcs.file.path, open = "rb")
  on.exit(close(con))
  ## seek the byte position where the DATA segment/stream starts
  seek(con, offsets['start'])
  ## read data stream: binary --> numeric --> matrix --> data.table
  dt <- data.table::as.data.table(
    matrix(
      data = readBin(
        con = con,
        what = "numeric",
        n = (diff(offsets) + 1)/(32/8),
        size = 32/8,
        signed = TRUE,
        endian = attributes(offsets)$endianness
      ),
      ncol = attributes(offsets)$par.n,
      byrow = TRUE,
      dimnames = list(NULL,(pars.N))
    )
  )
  ## return the data.table
  invisible(dt)
}

flowstate.from.file.path <- function(fcs.file.path){
  ## keyword-value pairs from TEXT segment; based on offsets from 'FCSoffsets(...)[["TEXT"]]'
  keywords <- readFCStext(fcs.file.path)
  ## create a flowstate (fs) S3 object ; class 'flowstate'
  fs <- flowstate(
    ## DATA segment; based on (updated) offsets from 'FCSoffsets(...)[["DATA"]]'
    data = readFCSdata(fcs.file.path),
    parameters = parameters.to.data.table(keywords, add.PROJ.identifier = TRUE),
    keywords = keywords.to.data.table(keywords, drop.primary = TRUE, drop.spill = TRUE),
    spill = spill.to.data.table(keywords)
  )
  ## any/all .fcs files should return:
  ## [['data']], [['parameters']], and [['keywords']]
  ## may not have spill (mass cytometry)
  if(fs[['spill']][,.N]==0){fs[['spill']] <- NULL}
  ##
  return(fs)
}

#' @title `flowstate`: read, parse, and store .fcs data
#' @description
#' [Flow Cytometry Standard](https://pmc.ncbi.nlm.nih.gov/articles/PMC2892967/) files are read, parsed and stored as objects (S3) of [class] `'flowstate'`.
#'
#' The individual segments (`TEXT` and `DATA`) of the FCS file are parsed and stored as follows:
#' * `[['data']]` -- a [data.table][data.table::data.table] containing raw/linear measurement values (scatter, MFI, Time, etc.).
#' * `[['parameters']]` -- a [data.table][data.table::data.table] containing instrument-specific parameters ('$PnN','$PnS', etc.).
#' * `[['keywords']]` -- a [data.table][data.table::data.table] containing instrument/sample-specific keyword-value pairs (metadata).
#' * `[['spill']]` -- a [data.table][data.table::data.table] containing (if present) spillover values.
#'
#' @details
#' Access the individual named list elements as follows (assuming `flowstate` object is named `fs`):
#' * `fs$data` or `fs[['data']]`
#' * `fs$parameters` or `fs[['parameters']]`
#' * `fs$keywords` or `fs[['keywords']]`
#' * `fs$spill` or `fs[['spill']]`
#'
#' Other `flowstate` functions operate on the entire object (`fs`) and access specific list elements as needed.
#'
#' @param fcs.file.paths Character string; path(s) returned from `list.files(...,full.names=T,pattern=".fcs")`.
#' @param colnames.type Character string; one of:
#' \itemize{
#'   \item `"N"` -- `[['data']]` columns are named by using only their respective $PN (name) keyword value.
#'   \item `"S"` -- `[['data']]` columns are named by using only their respective $PS (stain) keyword value.
#'   \item `"S_N"` -- `[['data']]` columns are named by combining $PS (stain) and $PN (name), separated by an underscore.
#' }
#' @param sample.id Keyword name -- default `NULL`; based on cytometer platform, the keyword name for `sample.id` will be automatically set and the respective keyword values will be added to `[['data']]` as a factored sample identifier. One of:
#' \itemize{
#'   \item Aurora (Cytek) -- `sample.id` = `'TUBENAME'`
#'   \item ID7000 (Sony)  -- `sample.id` = `'$CELLS'`
#'   \item Unspecified    -- `sample.id` = `'$FIL'`
#' }
#' @param concatenate Logical -- default `FALSE`; if `TRUE`, the list of flowstate objects will be combined into a single flowstate object.
#'
#' @returns For a single file: an object of class `flowstate`; for multiple files: a named list of `flowstate objects`; for concatenated files: an object of class `flowstate`
#' @seealso [flowstate.transform()]
#' @references Directly/heavily-inspired by:
#'
#' [flowCore][flowCore::flowFrame]:
#'
#' Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J, Jiang M, Finak G (2024). _flowCore: flowCore: Basic
#' structures for flow cytometry data_. doi:10.18129/B9.bioc.flowCore <https://doi.org/10.18129/B9.bioc.flowCore>, R package
#' version 2.22.0, <https://bioconductor.org/packages/flowCore>.
#'
#' [data.table][data.table::data.table]:
#'
#' Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T, Schwendinger B, Krylov I (2025).
#' data.table: Extension of 'data.frame'. R package version 1.17.99, https://r-datatable.com.
#'
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read a single .fcs file as a flowstate object
#' fs <- read.flowstate(
#'   fcs.file.paths[1],
#'   colnames.type = "S"
#' )
#' class(fs)
#' names(fs)
#'
#' #.fcs DATA segment as a data.table
#'   fs$data
#' #.fcs TEXT segment parsed and stored as three elements (data.tables):
#'   fs$parameters #instrument-specific measurement parameters
#'   fs$keywords #instrument/sample-specific metadata
#'   fs$spill #instrument/sample-specific spillover
#'
#' #read all .fcs files as a named list containing individual flowstate objects
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = FALSE
#' )
#' class(fs);class(fs[[1]])
#' names(fs);names(fs[[1]])
#'
#' fs[[1]]$keywords
#' fs[[1]]$data[, levels(sample.id)]
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#' class(fs)
#' names(fs)
#' fs$keywords
#' fs$data[, levels(sample.id)]
#'
read.flowstate<-function(
    fcs.file.paths,
    colnames.type = c("N", "S", "S_N"),
    sample.id = NULL,
    concatenate = FALSE
)
{
  ## add names to fcs.files.paths
  fcs.file.paths <- stats::setNames(fcs.file.paths, nm = basename(fcs.file.paths))
  ## create the flowstate object(s)
  fs <- lapply(fcs.file.paths, function(fcs.file.path){
    message(paste(basename(fcs.file.path), "-->", "flowstate"))
    flowstate.from.file.path(fcs.file.path)
  })
  ##
  colnames.type <- switch(
    match.arg(colnames.type),
    N = "N.alias",
    S = "S.alias",
    S_N = "S_N.alias"
  )
  ## update [['data']] names by reference using data.table::setnames
  ## update [['spill']] to match
  invisible(
    lapply(fs, function(.fs){
      if(all(names(.fs$data) != attributes(.fs$parameters)[['alias']][['N']])){
        stop("Naming conflict between [['data']] and [['parameters']]")
      }
      j <- attributes(.fs$parameters)[['alias']][[colnames.type]]
      data.table::setnames(x = .fs$data, old = names(.fs$data), new = j)
      if('spill' %in% names(.fs)){
        data.table::setnames(
          x = .fs$spill,
          new = j[match(names(.fs$spill), .fs$parameters[['N']])]
        )
      }
    })
  )
  ## get keyword identifier for use in defining 'sample.id' argument
  if(is.null(sample.id)){
    sample.id <- cytometer.identifier(fcs.file.paths)
  }
  ## add the identifier from [['keywords']] to [['keywords']] as sample.id (factor)
  ## add the factored identifier to [['data']]
  invisible(
    lapply(fs, function(.fs){
      if(!sample.id %in% names(.fs$keywords)){
        sample.id <- '$FIL'
      }
      data.table::set(
        x = .fs$keywords,
        j = 'sample.id',
        value = factor(sub(".fcs$", "", .fs$keywords[[sample.id]], ignore.case = T))
      )
      data.table::set(
        x = .fs$data,
        j = 'sample.id',
        value = .fs$keywords[, rep(sample.id, as.numeric(`$TOT`))]
      )
    })
  )
  ## return the object(s)
  if(length(fs) == 1){
    return(fs[[1]])
  }else{
    if(concatenate){
      fs <- concatenate.flowstate(fs)
    }else{
      return(fs)
    }
  }
}
