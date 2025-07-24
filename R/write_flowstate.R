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
write.flowstate<-function(flowstate.object,filename,endianness = c('little','big')){
  byteord<-switch(
    match.arg(endianness),
    little = "1,2,3,4",
    big = "4,3,2,1"
  )
  ##NULL identifier(s) in [['data']]
  ##identifiers will have a class of factor
  cols.null<-flowstate.object$data[,names(.SD),.SDcols = is.factor]
  flowstate.object$data[,(cols.null):=NULL]
  ##reverse any transformations
  flowstate.transform.inverse(flowstate.object)
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
    unlist(flowstate.object$keywords),
    spill.to.string(flowstate.object)
  )
  keyword.list<-as.list(keyword.string)
  keyword.list<-keyword.list[order(names(keyword.list))]
  ##update; flowstate.object-specific
  keyword.list[['$TOT']]<-as.character(flowstate.object$data[,.N])
  keyword.list[['$PAR']]<-as.character(ncol(flowstate.object$data))
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
  # for(i in c("BEGIN","END")){
  #   .pattern<-paste0(i,"DATA|0")
  #   regmatches(
  #     x = text.segment,
  #     m = gregexpr(.pattern,text.segment,fixed = TRUE)
  #   )<-paste(paste0(i,"DATA"),DATA.start,sep="|")
  # }
  # regmatches(
  #   x = text.segment,
  #   m = gregexpr("BEGINDATA|0",text.segment,fixed = TRUE)
  # )<-paste("BEGINDATA",DATA.start,sep = "|")
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
