FCSoffsets<-function(con){
  ##HEADER segment
  ##version identifier; TEXT segment; DATA segment; ANALYSIS segment
  ##6+4+8*2+8*2+8*2 = 58 bytes
  ##initialize a list for storing byte offsets
  offsets<-stats::setNames(
    vector(mode = "list",length = 4L),
    nm = c("version","TEXT","DATA","ANALYSIS")
  )
  ##first six bytes contain "FCS#.#"; version identifier
  version <- readChar(con, 6)
  if(version!="FCS3.1"){
    stop("Is this a .fcs file (version 3.1)?")
  }
  version <- substring(version, 4, nchar(version))
  offsets$version<-as.double(version)
  ##next four bytes contain space characters (ASCII 32);
  ##stop otherwise (not a legal/expected .FCS file)
  stopifnot(readChar(con,4) == "    ")
  ##byte offsets for TEXT, DATA, and ANALYSIS
  ##character strings; 8 bytes; right justified with space padding
  ##start and end positions
  for(i in c("TEXT","DATA","ANALYSIS")){
    offsets[[i]]<- c(
      start = as.integer(readChar(con,8)),
      end = as.integer(readChar(con,8))
    )
  }
  return(offsets)
}

readFCStext<-function(con,offsets){
  ##using the return of 'FCSoffsets(...)$TEXT'
  seek(con,offsets['start'])
  ##hex/ASCII encoded keyword-value pairs
  txt <- readBin(
    con = con,
    what = "raw",
    n = offsets['end'] + 1 - offsets['start']
  )
  ##convert
  txt<-rawToChar(txt)
  ##delimiter
  delimiter<-substr(txt, 1, 1)
  ##split the string; start at position 2
  kv<-strsplit(substr(txt,2,nchar(txt)),delimiter,fixed = TRUE)[[1]]#match exactly the delimiter
  ##keyword value pairs; keyword name (odd) and value (even)
  kv<-stats::setNames(
    kv[seq_along(kv) %%2 != 1],
    nm=kv[seq_along(kv) %%2 == 1]
  )
  ##as a list
  kv<-as.list(kv)
  return(kv)
}

parameters.to.data.table<-function(
    keywords,
    add.PROJ.identifier=TRUE,
    colnames.syntactically.valid=TRUE
)
{
  ##using the return of 'readFCStext(...)'
  ##number of parameters
  par.n<-as.numeric(keywords[['$PAR']])
  ##regex pattern to find '$P#' and 'P#'
  par.pattern<-"^\\$P\\d+|^P\\d+"
  parameter.names<-grep(
    pattern = par.pattern,
    x = names(keywords),
    value = TRUE
  )
  ##parameters; named vector
  pars<-trimws(unlist(keywords[parameter.names]))
  ##split parameters vector into 'types'; list
  par.types<-split(x=pars,f=factor(sub(par.pattern,"",names(pars))))
  ##create a data.table
  dt.parameters<-data.table::data.table(par=paste0('$P',seq(par.n)))
  ##fill the data.table by 'par.type'
  for(j in names(par.types)){
    par.vec<-par.types[[j]]
    par.i<-as.integer(gsub("\\D+","",names(par.vec)))
    data.table::set(dt.parameters,i=par.i,j=j,value = par.vec)
  }
  ##add '$PROJ' identifier; if not found, use '$DATE' instead
  if(add.PROJ.identifier){
    if(any(grepl("PROJ",names(keywords)))){
      proj<-grep("PROJ",names(keywords),value = T)
      dt.parameters[,PROJ:=as.factor(keywords[[proj]])]
    }else if(any(grepl("DATE",names(keywords)))){
      date<-grep("DATE",names(keywords),value = T)
      dt.parameters[,PROJ:=as.factor(keywords[[date]])]
    }
  }
  ##make N syntactically valid
  ##make S syntactically valid
  ##create alias columns
  if(colnames.syntactically.valid){
    dt.parameters[!grepl("[FS]SC|Time",N),N.alias:=gsub(" |-","",sub("-A$","",N))]
    dt.parameters[grepl("[FS]SC",N),N.alias:=sub("-","_",sub("SSC-B","SSCB",N))]
    dt.parameters[is.na(N.alias),N.alias := N]
    ##
    dt.parameters[,S.alias:=gsub(" |-","",S)]
    dt.parameters[is.na(S.alias),S.alias:=N.alias]
    ##
    dt.parameters[!is.na(S),S_N.alias := paste(S.alias,N.alias,sep="_")]
    dt.parameters[is.na(S_N.alias),S_N.alias := N.alias]
  }
  ##return the data.table
  dt.parameters[]
}

keywords.to.data.table<-function(keywords,drop.primary=TRUE,drop.spill=TRUE){
  ##using the return of 'readFCStext(...)'
  ##regex pattern to find '$P#' and 'P#'
  par.pattern<-"^\\$P\\d+|^P\\d+"
  keyword.names<-grep(
    pattern = par.pattern,
    x = names(keywords),
    value = TRUE,
    invert = T
  )
  ##keywords; named vector
  kw<-trimws(unlist(keywords[keyword.names]))
  ##create a data.table
  dt.keywords<-data.table::as.data.table(as.list(kw))
  ##drop required keywords
  if(drop.primary){
    dt.keywords[,(fcs.text.primary.required.keywords):=NULL]
  }
  ##drop spillover keyword-value pair (string)
  if(drop.spill){
    spill.name<-grep("spill",names(dt.keywords),ignore.case = T,value = T)
    dt.keywords[,(spill.name):=NULL]
  }
  ##add keywords related to originality/modification
  dt.keywords[,'$ORIGINALITY' := "DataModified"]
  dt.keywords[,'$LAST_MODIFIED' := toupper(format(Sys.time(), "%d-%b-%Y %H:%M:%OS2"))]
  dt.keywords[,'$LAST_MODIFIER' := Sys.getenv("USERNAME")]

  ##return the data.table
  dt.keywords[]
}

spill.to.data.table<-function(keywords,add.PROJ.identifier=TRUE){
  ##using the return of 'readFCStext(...)'
  ##spill/spillover keyword name
  spill.name<-grep('spill',names(keywords),ignore.case = TRUE,value = T)
  ##parse the spillover string
  if(length(spill.name)==1){
    spillover.string<-keywords[[spill.name]]
    spill.split<-unlist(strsplit(spillover.string,","))
    N.cols<-as.numeric(spill.split[1])
    col.names <- spill.split[2:(N.cols + 1)]
    vals<-as.numeric(spill.split[(N.cols+2):length(spill.split)])
    dt.spill<-data.table::as.data.table(
      matrix(
        data = vals,
        ncol = N.cols,
        byrow = TRUE,
        dimnames = list(NULL,col.names)
      )
    )
    rownames(dt.spill)<-names(dt.spill)
    if(add.PROJ.identifier){
      proj<-grep("PROJ",names(keywords),value = T)
      data.table::setattr(dt.spill,'PROJ',keywords[[proj]])
      # attr(spill.mat,'PROJ')<-keywords[[proj]]
      # return(spill.mat)
    }
    ##
    data.table::setattr(dt.spill,name = "applied",value = FALSE)
    ##
    return(dt.spill)
  }
}

offsets.data.update<-function(offsets,keywords){
  ##update offsets$DATA' with '$BEGINDATA' and '$ENDDATA' keyword values
  ##add 'endianness' based on keyword '$BYTEORD'
  ##add 'par.n' based on keyword '$PAR'
  offsets<-replace(
    offsets,
    values=as.numeric(unlist(keywords[c('$BEGINDATA','$ENDDATA')]))
  )
  endian.string<-c("4,3,2,1"='big',"1,2,3,4"='little')
  if(!keywords$`$BYTEORD` %in% names(endian.string)){
    stop("Keyword '$BYTEORD' has unexpected endianness...data segment cannot be appropriately parsed.")
  }else{
    attr(offsets,"endianness") <- endian.string[[keywords[['$BYTEORD']]]]
    attr(offsets,"par.n") <- as.numeric(keywords[['$PAR']])
  }
  ##
  return(offsets)
}

readFCSdata<-function(con,offsets){
  ##using the return of 'FCSoffsets(...)$DATA' --> updated with 'update.offsets.data(...)'
  ##DATA stream start
  seek(con,offsets['start'])
  ##read data stream --> matrix --> data.table
  dt<-data.table::as.data.table(
    matrix(
      data = readBin(
        con = con,
        what = "numeric",
        n = (offsets[['end']] + 1 - offsets[['start']])/(32/8),
        size = 32/8,
        signed = TRUE,
        endian = attributes(offsets)$endianness
      ),
      ncol = attributes(offsets)$par.n,
      byrow = TRUE
    )
  )
  ##return the data.table
  dt[]
}

flowstate.object<-function(fcs.file){
  ##open a connection to the file (.fcs); read binary mode
  con <- file(fcs.file, open = "rb")
  on.exit(close(con))
  ##get byte offsets for reading FCS version, TEXT, DATA, and ANALYSIS segments
  offsets<-FCSoffsets(con)
  ##keyword-value pairs from TEXT segment
  keywords<-readFCStext(con,offsets$TEXT)
  ##update offsets$DATA
  offsets$DATA<-offsets.data.update(offsets$DATA,keywords)
  ##create a flowstate.object (fs); class 'flowstate'
  fs<-structure(
    list(
      data = readFCSdata(con,offsets$DATA),
      parameters = parameters.to.data.table(keywords,add.PROJ.identifier = TRUE),
      keywords = keywords.to.data.table(keywords,drop.primary = TRUE,drop.spill = TRUE),
      spill = spill.to.data.table(keywords),
      meta = NULL
    ),
    class = "flowstate"
  )
  ##
  return(fs)
}

#' @title flowstate: read, process, and store .fcs data
#'
#' @param fcs.file.paths Character string; path(s) returned from `list.files(...,full.names=T,pattern=".fcs")`.
#' @param colnames.type Character string; one of:
#' \itemize{
#'   \item `"S_N"` -- `[['data']]` columns are named by combining $PS (stain) and $PN (name), separated by an underscore.
#'   \item `"S"` -- `[['data']]` columns are named by using only their respective $PS (stain) keyword value.
#'   \item `"N"` -- `[['data']]` columns are named by using only their respective $PN (name) keyword value.
#' }
#' @param cofactor Numeric; default `5000`. Any/all parameters with a `$PnTYPE` of 'Raw_Fluorescence' or 'Unmixed_Fluorescence' will be transformed using \link{asinh} and the defined cofactor value (`asinh(x/cofactor)`).
#' @param sample.id Character string; the keyword label defined through `sample.id` (default `TUBENAME`) will be used to add respective keyword values from `[['keywords']]` as an identifier to `[['data']]`.
#'
#' @returns An object of class flowstate.
#' @export
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = ".fcs")
#'
#' #read a single .fcs file as a flowstate object
#' fs <- read.flowstate(
#'   fcs.file.paths[1],
#'   colnames.type="S",
#'   cofactor = 5000
#' )
#' class(fs)
#'
#' #.fcs DATA segment as a data.table
#'   fs$data
#' #.fcs TEXT segment parsed and stored as three elements (data.tables)
#'   fs$parameters #instrument-specific measurement parameters
#'   fs$keywords #instrument/sample-specific metadata
#'   fs$spill #instrument/sample-specific spillover
#'
#' #plot some data
#' no.fill.legend <- ggplot2::guides(fill = 'none')
#' .title1 <- paste("Batch:",fs$keywords[,`$PROJ`])
#' .title2 <- paste(
#'   "Instrument Serial#:",
#'   fs$keywords[,paste(.(`$CYT`,`$CYTSN`),collapse = " ")]
#' )
#' .title <- paste(.title1,.title2,sep = "\n")
#' .title <- ggplot2::labs(title = .title, subtitle = fs$keywords[,TUBENAME])
#'
#' plot(fs,CD3,Viability) + no.fill.legend + .title
#' plot(fs,FSC_A,SSC_A) + no.fill.legend + .title
#'
read.flowstate<-function(
    fcs.file.paths,
    colnames.type=c("S_N","S","N"),
    cofactor=5000,
    sample.id='TUBENAME'
)
{
  ##add names to fcs.files.paths
  fcs.file.paths<-stats::setNames(fcs.file.paths,nm=basename(fcs.file.paths))
  ##create the object(s)
  fs<-lapply(stats::setNames(fcs.file.paths,nm=basename(fcs.file.paths)),flowstate.object)
  ##
  colnames.type<-switch(
    match.arg(colnames.type),
    S_N = "S_N.alias",
    S = "S.alias",
    N = "N.alias"
  )
  ##update [['data']] by reference using data.table::setnames
  invisible(
    lapply(fs,function(fs.obj){
      data.table::setnames(fs.obj$data,new = fs.obj$parameters[[colnames.type]])
    })
  )
  ##update [['spill']] to match
  invisible(
    lapply(fs,function(fs.obj){
      data.table::setnames(
        x = fs.obj$spill,
        new = fs.obj$parameters[[colnames.type]][match(names(fs.obj$spill),fs.obj$parameters[['N']])]
      )
    })
  )
  ##transform [['data']]
  if(!is.null(cofactor)){
    invisible(
      lapply(fs,function(fs.obj){
        ##Cytek Aurora; spectral
        cols.transform<-fs.obj$parameters[grep("fluorescence",TYPE,ignore.case = T)][[colnames.type]]
        ##
        for(j in cols.transform){
          data.table::set(
            x = fs.obj$data,
            j = j,
            value = asinh(fs.obj$data[[j]]/cofactor)
          )
        }
        ##
        fs.obj$parameters[
          i = fs.obj$parameters[[colnames.type]] %in% cols.transform,
          j = c('transform','cofactor'):=list('asinh',cofactor)
        ]
      })
    )
  }
  ##add an identifier to [['data']]
  invisible(
    lapply(fs,function(fs.obj){
      if(!sample.id %in% names(fs.obj$keywords)){
        sample.id<-'$FIL'
      }
      data.table::set(
        x = fs.obj$data,
        j = 'sample.id',
        value=as.factor(
          fs.obj$keywords[
            ,
            rep(sub(".fcs","",j),as.numeric(`$TOT`)),
            env = list(j = sample.id)
          ]
        )
      )
    })
  )
  ##return the object(s)
  if(length(fs)==1){
    return(fs[[1]])
  }else{
    return(fs)
  }
}
