parms.update<-function(flowstate){
  j.match <- j.match.parameters.to.data(flowstate)
  parms.update <- names(flowstate$data)[
    !names(flowstate$data) %in% flowstate$parameters[[j.match]]]
  parms.update <- parms.update[!parms.update %in% 'sample.id']
}
data.cols.to.numeric <- function(flowstate){
  parms.update.convert<-names(which(flowstate$data[
    ,
    sapply(.SD,function(j){class(j)!='numeric'}),
    .SDcols = parms.update(flowstate)]))
  ##convert to numeric in [['data']]
  if(length(parms.update.convert)>=1){
    message(paste("Converting",paste0(parms.update.convert,collapse = ", "), "to numeric"))
    message("Updating [['data']] by reference using data.table::set()")
    for(j in parms.update.convert){
      data.table::set(
        x = flowstate$data,
        j = j,
        value = as.numeric(flowstate$data[[j]])#returns factor levels
      )
    }
  }else{
    message("No columns identified for conversion to numeric.")
  }
}
parameters.dt.update<-function(flowstate){
  .parms.update<-parms.update(flowstate)
  ##create a [['parameters']] data.table using any/all parms in .parms.update
  dt.parms.update<-data.table::data.table(
    par = paste0('$P',seq(length(.parms.update))+flowstate$parameters[,.N]),
    B = "32",
    DISPLAY = "LIN",
    E = "0,0",
    N = .parms.update,
    R = flowstate$data[
      ,
      sapply(.SD,function(j){as.character(max(ceiling(as.numeric(j))))}),
      .SDcols = .parms.update],
    TYPE = "Derived_Numeric"
  )
  ##updated [['parameters']]; needs reassignment
  dt.parms.update <- rbind(
    flowstate$parameters,
    dt.parms.update,
    fill = TRUE
  )
}
keywords.to.character<-function(flowstate){
  keywords.convert<-names(which(flowstate$keywords[
    ,
    sapply(.SD,function(j){class(j)!='character'})]))
  ##convert to character in [['keywords']]
  if(length(keywords.convert)>=1){
    message(paste("Converting",paste0(keywords.convert,collapse = ", "), "to character."))
    message("Updating [['keywords']] by reference using data.table::set().")
    for(j in keywords.convert){
      data.table::set(
        x = flowstate$keywords,
        j = j,
        value = as.character(flowstate$keywords[[j]])
      )
    }
  }else{
    message("No columns identified for conversion to character.")
  }
}
##a few checks with some verbose output;
##eventually update with automatic execution of the individually named functions within the check
write.flowstate.check <- function(flowstate,verbose=TRUE){
  ##
  j.match <- j.match.parameters.to.data(flowstate)
  ##
  if('transform' %in% names(flowstate$parameters)){
    res <- flowstate$parameters[!is.na(transform)][[j.match]]
    if(length(res)>0){
      if(verbose){
        message("Inverse transform the following columns in [['data']] using flowstate:::flowstate.transform.inverse():")
        print(res)
      }
    }
  }
  ##
  res<-names(which(flowstate$data[,sapply(.SD,class)!="numeric",.SDcols = !'sample.id']))
  if(length(res)>0){
    if(verbose){
      message("Convert the following columns in [['data']] to numeric using flowstate:::data.cols.to.numeric():")
      print(res)
    }
  }
  ##this check needs fixing! Returns NAs after the update
  res <- names(flowstate$data)[!names(flowstate$data) %in% flowstate$parameters[[j.match]]]
  res <- grep("sample.id",res,value = T,invert = T)
  if(length(res)>0){
    if(verbose){
      message("Re-assign [['parameters']] <- flowstate:::parameters.dt.update() to include:")
      print(res)
    }
  }
  ##
  res <- names(which(flowstate$keywords[,sapply(.SD,class)=="factor"]))
  if(length(res)>0){
    if(verbose){
      message("Convert the following columns in [['keywords']] to character using flowstate:::keywords.to.character():")
      print(res)
    }
  }
  ##
  res <- attr(flowstate$spill,'applied')
  if(res){
    if(verbose){
      message("Decompensate using flowstate::spillover.apply(...,decompensate = TRUE)")
      print(res)
    }
  }
}
##flowstate parameters (data.table) to vector (string); function
parameters.to.string<-function(flowstate){
  cols.parameters<-names(flowstate$parameters)[
    names(flowstate$parameters) %in% flowstate.parameter.keywords
  ]
  parms<-lapply(cols.parameters,function(p,par.n=flowstate$parameters[,.N]){
    vec<-stats::setNames(
      as.character(flowstate$parameters[[p]]),
      nm=paste0(ifelse(p=='DISPLAY',"P","$P"),seq(par.n),p)
    )
    return(vec)
  })
  parms<-unlist(parms)
  parms<-parms[order(names(parms))]
  parms<-parms[!is.na(parms)]
  ##
  return(parms)
}
##flowstate spill (data.table) to vector (string); function
spill.to.string<-function(fs.obj){
  spill<-data.table::copy(fs.obj$spill)
  col.match<-fs.obj$parameters[,sapply(.SD,function(j){all(names(spill) %in% j)})]
  col.match<-names(which(col.match))
  if(length(col.match)>1) col.match<-col.match[1]
  names.match<-fs.obj$parameters[['N']][match(names(spill),fs.obj$parameters[[col.match]])]
  data.table::setnames(spill,new = names.match)
  ##
  spill.string<-paste0(
    c(
      spill[,.N],
      paste0(names(spill),collapse = ","),
      paste0(c(t(as.matrix(spill))),collapse = ",")
    ),
    collapse = ","
  )
  names(spill.string)<-'$SPILLOVER'
  ##
  return(spill.string)
}
#' @title Write a `flowstate` as a Flow Cytometry Standard file (FCS 3.1)
#' @description
#' A `flowstate` -- following [Flow Cytometry Standard](https://pmc.ncbi.nlm.nih.gov/articles/PMC2892967/) conventions -- is converted to FCS 3.1 and written to disk. The conversion is as follows:
#' * `[['keywords']]`, `[['parameters']]`, and `[['spill']]` (if present) are converted to character string using a delimiter `"|"` to separate keyword-value pairs.
#'   * this string forms the `TEXT` segment of the FCS file.
#' * `[['data']]` is written as binary.
#'   * this forms the `DATA` segment of the FCS file.
#' * The `HEADER` segment (string) -- containing the required offsets -- is derived from `TEXT` and `DATA`.
#'
#' The FCS 3.1 compliant file is then written to disk.
#'
#' *N.B.: During testing, FCS 3.1 files written to disk by [write.flowstate] were fully compatible with [read.flowstate] and [flowCore::read.FCS]; issues were only with proprietary software (i.e., FlowJo). If proprietary software is to be used with FCS 3.1 files generated by [write.flowstate], caution is advised.*
#'
#' @param flowstate A `flowstate`
#' @param new.fil Character string -- default `NULL`; if defined, the keyword `'$FIL'` will be updated with the provided value.
#' @param add.fil.mod Logical -- default `TRUE`; the suffix `'flowstateMOD'`(string) will be appended to the value of `'$FIL'`.
#' @param file.dir Character vector -- a [base][file.path]; a file directory for saving the output FCS 3.1 file.
#' @param endianness Character string -- default `"little"`; see [endian][base::writeBin].
#'
#' @returns The FCS 3.1 file is written to disk and a summary message is displayed.
#' @keywords internal
#'
#' @examples
#'
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
#'
#' #read .fcs files as a flowstate object; concatenate
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type = "S",
#'   concatenate = TRUE
#' )
#'
#' #remove saturating events
#' select_nonsaturating(fs)
#'
#' #subset to retain only non-saturating events
#' fs <- subset(fs, select.nonsaturating)
#' fs$data[, select.nonsaturating := NULL]
#'
#' #transform
#' flowstate.transform(
#'   fs,
#'   j = c('CD3','CD4','CD8'),
#'   transform.func = "asinh",
#'   cofactor = 5000
#' )
#' #
write.flowstate <- function(
    flowstate,
    new.fil = NULL,
    add.fil.mod = TRUE,
    file.dir,
    endianness = c('little', 'big')
)
{
  ## PREPARE flowstate START
  ## inverse transformation of [['data']] columns
    flowstate.transform.inverse(flowstate)

  # ## PREPARE flowstate START
  # j.match <- j.match.parameters.to.data(flowstate)

  # ## inverse transformation of [['data']] columns
  # if('transform' %in% names(flowstate$parameters)){
  #   res <- flowstate$parameters[,any(!is.na(transform))]
  #   if(res){
  #     flowstate.transform.inverse(flowstate)
  #   }
  # }

  ##convert non-numeric (factored) [['data']] columns
  res <- flowstate$data[
    ,
    any(sapply(.SD,class)!="numeric"),
    .SDcols = !'sample.id'
  ]
  if(res){
    suppressMessages(
      data.cols.to.numeric(flowstate)
    )
  }
  ##convert non-character (factored/numeric) [['keyword']] columns
  res <- any(flowstate$keywords[,sapply(.SD,class)!="character"])
  if(res){
    suppressMessages(
      keywords.to.character(flowstate)
    )
  }
  ##update [['parameters']] to match [['data']]; derived columns
  res <-any(
    !flowstate$data[,names(.SD),.SDcols = !'sample.id'] %in%
      flowstate$parameters[[j.match]]
  )
  if(res){
    ##needs re-assignment
    flowstate$parameters <- parameters.dt.update(flowstate)
  }
  ##PREPARE flowstate END

  ##
  byteord<-switch(
    match.arg(endianness),
    little = "1,2,3,4",
    big = "4,3,2,1"
  )
  ##NULL identifier(s) in [['data']]
  ##identifiers will have a class of factor
  cols.null<-flowstate$data[,names(.SD),.SDcols = is.factor]
  flowstate$data[,(cols.null):=NULL]
  ##required FCS keywords
  keywords.required<-fcs.text.primary.required.keywords
  keywords.required<-stats::setNames(nm=keywords.required,rep("0",length(keywords.required)))
  keywords.required[['$BYTEORD']]<-byteord
  keywords.required[['$DATATYPE']]<-'F'
  keywords.required[['$MODE']]<-"L"
  ##vector (string) of keywords
  keyword.string<-c(
    keywords.required,
    parameters.to.string(flowstate),
    unlist(flowstate$keywords)
    # flowstate:::spill.to.string(flowstate)
  )
  if('spill' %in% names(flowstate)){
    keyword.string<-c(
      keyword.string,
      spill.to.string(flowstate)
    )
  }
  keyword.list<-as.list(keyword.string)
  keyword.list<-keyword.list[order(names(keyword.list))]
  ##update; flowstate-specific
  keyword.list[['$TOT']]<-as.character(flowstate$data[,.N])
  keyword.list[['$PAR']]<-as.character(ncol(flowstate$data))
  ##update '$FIL'
  if(!is.null(new.fil)){
    if(add.fil.mod){
      new.fil <- sprintf("%s_flowstateMOD",new.fil)
    }
    if(!grepl(".fcs$",new.fil)){
      new.fil <- paste0(new.fil,".fcs")
    }
    keyword.list[['$FIL']]<-new.fil
  }
  ##TEXT segment
  text.segment<-paste0(
    "|",
    paste0(names(keyword.list),"|",keyword.list,"|",collapse = "")
  )
  TEXT.start<-58
  TEXT.end<-nchar(text.segment, "bytes") + TEXT.start - 1

  data.stream.bytes <- flowstate$data[,.N] * ncol(flowstate$data) * 4

  ##from flowCore::write.FCS()
  ##in 'text.segment': $BEGINDATA = "0" and $ENDDATA = "0"
  ##combined, nchar = 2
  kw.len.old <- 2
  repeat {
    DATA.start <- TEXT.end + 1
    DATA.end <- DATA.start + data.stream.bytes - 1
    kw.len.new <- nchar(DATA.start) + nchar(DATA.end)
    if (kw.len.new > kw.len.old) {
      TEXT.end <- TEXT.end + kw.len.new - kw.len.old
      kw.len.old <- kw.len.new
    }
    else break
  }

  ##update 'text.segment'
  regmatches(text.segment,gregexpr("BEGINDATA|0",text.segment,fixed = TRUE))<-paste("BEGINDATA",DATA.start,sep = "|")
  regmatches(text.segment,gregexpr("ENDDATA|0",text.segment,fixed = TRUE))<-paste("ENDDATA",DATA.end,sep = "|")
  ##test
  if(nchar(text.segment, "bytes") + TEXT.start - 1 != DATA.start-1){
    stop("Issue with updating 'text.segment' with '$BEGINDATA' and/or '$ENDDATA' values.")
  }
  ##initialize a list for storing byte offsets
  offsets<-stats::setNames(
    vector(mode = "list",length = 3L),
    nm = c("TEXT","DATA","ANALYSIS")
  )
  offsets$TEXT<-c(TEXT.start,TEXT.end)
  offsets$DATA<-c(DATA.start,DATA.end)
  offsets$ANALYSIS<-c(0,0)
  ##prepare the HEADER segment
  header.segment<-paste0(
    sapply(offsets,function(i){
      if(any(nchar(i)>8)){
        i<-c(0,0)
      }
      paste0(sprintf("%8s",i),collapse = "")
    }),
    collapse = ""
  )
  version<-sprintf("%-10s","FCS3.1")
  header.segment<-paste0(version,header.segment)
  if(nchar(header.segment) != 58){
    stop("HEADER segment does not match expected byte length/nchar of 58")
  }
  ##
  if(!dir.exists(file.dir)){dir.create(file.dir,recursive = T)}
  filename<-file.path(file.dir,keyword.list[['$FIL']])
  con <- file(filename, open = "wb")
  on.exit(close(con))
  seek(con,0)
  writeChar(header.segment, con, eos = NULL)
  writeChar(text.segment, con, eos = NULL)

  data.segment<-methods::as(t(flowstate$data),"numeric")
  writeBin(
    object = data.segment,
    con,
    size = 4,
    endian = match.arg(endianness)
  )
  writeChar("00000000", con, eos = NULL)
  ##
  obj <- deparse(substitute(flowstate))
  message(sprintf("%s --> %s",obj,filename))
}
