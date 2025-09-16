parms.update<-function(flowstate.object){
  j.match <- j.match.parameters.to.data(flowstate.object)
  parms.update <- names(flowstate.object$data)[
    !names(flowstate.object$data) %in% flowstate.object$parameters[[j.match]]]
  parms.update <- parms.update[!parms.update %in% 'sample.id']
}
data.cols.to.numeric <- function(flowstate.object){
  parms.update.convert<-names(which(flowstate.object$data[
    ,
    sapply(.SD,function(j){class(j)!='numeric'}),
    .SDcols = parms.update(flowstate.object)]))
  ##convert to numeric in [['data']]
  if(length(parms.update.convert)>=1){
    message(paste("Converting",paste0(parms.update.convert,collapse = ", "), "to numeric"))
    message("Updating [['data']] by reference using data.table::set()")
    for(j in parms.update.convert){
      data.table::set(
        x = flowstate.object$data,
        j = j,
        value = as.numeric(flowstate.object$data[[j]])#returns factor levels
      )
    }
  }else{
    message("No columns identified for conversion to numeric.")
  }
}
parameters.dt.update<-function(flowstate.object){
  .parms.update<-parms.update(flowstate.object)
  ##create a [['parameters']] data.table using any/all parms in .parms.update
  dt.parms.update<-data.table::data.table(
    par = paste0('$P',seq(length(.parms.update))+flowstate.object$parameters[,.N]),
    B = "32",
    DISPLAY = "LIN",
    E = "0,0",
    N = .parms.update,
    R = flowstate.object$data[
      ,
      sapply(.SD,function(j){as.character(max(ceiling(as.numeric(j))))}),
      .SDcols = .parms.update],
    TYPE = "Derived_Numeric"
  )
  ##updated [['parameters']]; needs reassignment
  dt.parms.update <- rbind(
    flowstate.object$parameters,
    dt.parms.update,
    fill = TRUE
  )
}
keywords.to.character<-function(flowstate.object){
  keywords.convert<-names(which(flowstate.object$keywords[
    ,
    sapply(.SD,function(j){class(j)!='character'})]))
  ##convert to character in [['keywords']]
  if(length(keywords.convert)>=1){
    message(paste("Converting",paste0(keywords.convert,collapse = ", "), "to character."))
    message("Updating [['keywords']] by reference using data.table::set().")
    for(j in keywords.convert){
      data.table::set(
        x = flowstate.object$keywords,
        j = j,
        value = as.character(flowstate.object$keywords[[j]])
      )
    }
  }else{
    message("No columns identified for conversion to character.")
  }
}
##a few checks with some verbose output;
##eventually update with automatic execution of the individually named functions within the check
write.flowstate.check <- function(flowstate.object,verbose=TRUE){
  ##
  j.match <- j.match.parameters.to.data(flowstate.object)
  ##
  res <- flowstate.object$parameters[!is.na(transform)][[j.match]]
  if(length(res)>0){
    if(verbose){
      message("Inverse transform the following columns in [['data']] using flowstate:::flowstate.transform.inverse():")
      print(res)
    }
  }
  ##
  res<-names(which(flowstate.object$data[,sapply(.SD,class)!="numeric",.SDcols = !'sample.id']))
  if(length(res)>0){
    if(verbose){
      message("Convert the following columns in [['data']] to numeric using flowstate:::data.cols.to.numeric():")
      print(res)
    }
  }
  ##this check needs fixing! Returns NAs after the update
  res <- names(flowstate.object$data)[!names(flowstate.object$data) %in% flowstate.object$parameters[[j.match]]]
  res <- grep("sample.id",res,value = T,invert = T)
  if(length(res)>0){
    if(verbose){
      message("Re-assign [['parameters']] <- flowstate:::parameters.dt.update() to include:")
      print(res)
    }
  }
  ##
  res <- names(which(flowstate.object$keywords[,sapply(.SD,class)=="factor"]))
  if(length(res)>0){
    if(verbose){
      message("Convert the following columns in [['keywords']] to character using flowstate:::keywords.to.character():")
      print(res)
    }
  }
  ##
  res <- attr(flowstate.object$spill,'applied')
  if(res){
    if(verbose){
      message("Decompensate using flowstate::spillover.apply(...,decompensate = TRUE)")
      print(res)
    }
  }
}
##flowstate parameters (data.table) to vector (string); function
parameters.to.string<-function(fs.obj){
  cols.parameters<-names(fs.obj$parameters)[
    names(fs.obj$parameters) %in% flowstate.parameter.keywords
  ]
  parms<-lapply(cols.parameters,function(p,par.n=fs.obj$parameters[,.N]){
    vec<-stats::setNames(
      fs.obj$parameters[[p]],
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
write.flowstate <- function(
    flowstate.object,
    new.fil=NULL,
    add.fil.mod=TRUE,
    file.dir,
    endianness = c('little','big')
)
{
  ##PREPARE flowstate.object START
  j.match <- j.match.parameters.to.data(flowstate.object)
  ##inverse transformation of [['data']] columns
  res <- flowstate.object$parameters[,any(!is.na(transform))]
  if(res){
    flowstate.transform.inverse(flowstate.object)
  }
  ##convert non-numeric (factored) [['data']] columns
  res <- flowstate.object$data[
    ,
    any(sapply(.SD,class)!="numeric"),
    .SDcols = !'sample.id'
  ]
  if(res){
    suppressMessages(
      data.cols.to.numeric(flowstate.object)
    )
  }
  ##convert non-character (factored/numeric) [['keyword']] columns
  res <- any(flowstate.object$keywords[,sapply(.SD,class)!="character"])
  if(res){
    suppressMessages(
      keywords.to.character(flowstate.object)
    )
  }
  ##update [['parameters']] to match [['data']]; derived columns
  res <-any(
    !flowstate.object$data[,names(.SD),.SDcols = !'sample.id'] %in%
      flowstate.object$parameters[[j.match]]
  )
  if(res){
    ##needs re-assignment
    flowstate.object$parameters <- parameters.dt.update(flowstate.object)
  }
  ##PREPARE flowstate.object END

  ##
  byteord<-switch(
    match.arg(endianness),
    little = "1,2,3,4",
    big = "4,3,2,1"
  )
  ##NULL identifier(s) in [['data']]
  ##identifiers will have a class of factor
  cols.null<-flowstate.object$data[,names(.SD),.SDcols = is.factor]
  flowstate.object$data[,(cols.null):=NULL]
  ##required FCS keywords
  keywords.required<-fcs.text.primary.required.keywords
  keywords.required<-stats::setNames(nm=keywords.required,rep("0",length(keywords.required)))
  keywords.required[['$BYTEORD']]<-byteord
  keywords.required[['$DATATYPE']]<-'F'
  keywords.required[['$MODE']]<-"L"
  ##vector (string) of keywords
  keyword.string<-c(
    keywords.required,
    parameters.to.string(flowstate.object),
    unlist(flowstate.object$keywords)
    # flowstate:::spill.to.string(flowstate.object)
  )
  if('spill' %in% names(flowstate.object)){
    keyword.string<-c(
      keyword.string,
      spill.to.string(flowstate.object)
    )
  }
  keyword.list<-as.list(keyword.string)
  keyword.list<-keyword.list[order(names(keyword.list))]
  ##update; flowstate.object-specific
  keyword.list[['$TOT']]<-as.character(flowstate.object$data[,.N])
  keyword.list[['$PAR']]<-as.character(ncol(flowstate.object$data))
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

  data.stream.bytes <- flowstate.object$data[,.N] * ncol(flowstate.object$data) * 4

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

  data.segment<-methods::as(t(flowstate.object$data),"numeric")
  writeBin(
    object = data.segment,
    con,
    size = 4,
    endian = match.arg(endianness)
  )
  writeChar("00000000", con, eos = NULL)
}
