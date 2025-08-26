FCSoffsets<-function(fcs.file.path){
  if(grepl(".fcs$",fcs.file.path) & file.exists(fcs.file.path)){
    ##open a connection to the file (.fcs); read binary mode
    con <- file(fcs.file.path, open = "rb")
    on.exit(close(con))
    ##first six bytes contain "FCS#.#"; version identifier
    version <- readChar(con, 6)
    ##version test
    if(!version %in% sprintf("FCS3.%d",0:1)){#sprintf("FCS3.%d",0:2)
      stop("Is this a .fcs file (version 3.0 or 3.1)? flowstate has only been tested using FCS 3.0/3.1 files.")
    }
  }else{
    stop("Need a valid '.fcs' file path/file.")
  }
  ##HEADER segment
  ##version identifier; TEXT segment; DATA segment; ANALYSIS segment
  ##6+4+8*2+8*2+8*2 = 58 bytes
  ##initialize a list for storing byte offsets
  offsets<-stats::setNames(
    vector(mode = "list",length = 4L),
    nm = c("version","TEXT","DATA","ANALYSIS")
  )
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

readFCStext<-function(fcs.file.path,con=NULL,offsets=NULL){
  ##use the return of 'FCSoffsets(...)[["TEXT"]]'
  offsets<-FCSoffsets(fcs.file.path)[["TEXT"]]
  ##open a connection to the file (.fcs); read binary mode
  con <- file(fcs.file.path, open = "rb")
  on.exit(close(con))
  ##seek the byte position where the TEXT segment starts
  seek(con, where = offsets['start'])
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

add.identifier.proj<-function(keywords){
  if(any(grepl("^\\$PROJ$",names(keywords)))){
    proj<-grep("^\\$PROJ$",names(keywords),value = T)
  }else if(any(grepl("$DATE",names(keywords),fixed = T))){
    proj<-grep("$DATE",names(keywords),value = T,fixed = T)
  }
  ##return the identifier
  keywords[[proj]]
}

parameters.to.data.table<-function(
    keywords,
    add.PROJ.identifier=TRUE,
    colnames.syntactically.valid=TRUE,
    S.func=NULL#function(j){strsplit(j," ")[[1]][1]}
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
    dt.parameters[,PROJ:=as.factor(add.identifier.proj(keywords))]
  }
  ##make N syntactically valid
  ##make S syntactically valid
  ##create alias columns
  if(colnames.syntactically.valid){
    dt.parameters[!grepl("[FS]SC|Time",N),N.alias:=gsub(" |-","",sub("-A$","",N))]
    dt.parameters[grepl("[FS]SC",N),N.alias:=sub("-","_",sub("SSC-B","SSCB",N))]
    dt.parameters[is.na(N.alias),N.alias := N]
    ##
    if("S" %in% names(dt.parameters)){
      if(is.null(S.func)){
        dt.parameters[,S.alias:=gsub(" |-","",S)]
      }else{
        dt.parameters[,S.alias := sapply(S,S.func)]
      }
      dt.parameters[is.na(S.alias),S.alias:=N.alias]
      ##
      dt.parameters[!is.na(S),S_N.alias := paste(S.alias,N.alias,sep="_")]
      dt.parameters[is.na(S_N.alias),S_N.alias := N.alias]
    }
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
    ##add '$PROJ' identifier; if not found, use '$DATE' instead
    if(add.PROJ.identifier){
      data.table::setattr(dt.spill,'PROJ',add.identifier.proj(keywords))
    }
    ##
    data.table::setattr(dt.spill,name = "applied",value = FALSE)
    ##
    return(dt.spill)
  }
}

offsets.data.update<-function(offsets,keywords){
  ##update offsets$DATA with '$BEGINDATA' and '$ENDDATA' keyword values
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

readFCSdata<-function(fcs.file.path,con=NULL,offsets=NULL){
  ##use the return of 'FCSoffsets(...)[["DATA"]]'
  offsets<-FCSoffsets(fcs.file.path)[["DATA"]]
  ##use the return of 'readFCStext(...)'
  keywords<-readFCStext(fcs.file.path)
  ##update 'offsets'
  offsets<-offsets.data.update(offsets,keywords)
  ##open a connection to the file (.fcs); read binary mode
  con <- file(fcs.file.path, open = "rb")
  on.exit(close(con))
  ##seek the byte position where the DATA segment/stream starts
  seek(con,offsets['start'])
  ##read data stream: binary --> numeric --> matrix --> data.table
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

flowstate.from.file.path<-function(fcs.file.path,S.func=NULL){
  ##keyword-value pairs from TEXT segment; based on offsets from 'FCSoffsets(...)[["TEXT"]]'
  keywords<-readFCStext(fcs.file.path)
  ##create a flowstate (fs) S3 object ; class 'flowstate'
  fs<-flowstate(
    ##DATA segment; based on (updated) offsets from 'FCSoffsets(...)[["DATA"]]'
    data = readFCSdata(fcs.file.path),
    parameters = parameters.to.data.table(keywords,S.func = S.func,add.PROJ.identifier = TRUE),
    keywords = keywords.to.data.table(keywords,drop.primary = TRUE,drop.spill = TRUE),
    spill = spill.to.data.table(keywords)
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
#' @param S.func a function; default `NULL`. If a function is supplied, it will be used to modify/split `"S"`; e.g. `function(j){strsplit(j," ")[[1]][1]}` will be applied to `"S"` to return the first split element ("CD4 PE" --> "CD4").
#' @param cofactor Numeric; default `5000`. Any/all parameters with a `$PnTYPE` of 'Raw_Fluorescence' or 'Unmixed_Fluorescence' will be transformed using \link{asinh} and the defined cofactor value (`asinh(x/cofactor)`).
#' @param sample.id Character string; the keyword label defined through `sample.id` (default `TUBENAME`) will be used to add respective keyword values from `[['keywords']]` as an identifier to `[['data']]`.
#' @param concatenate Logical; if `TRUE`, the list of flowstate objects will be combined into a single flowstate object.
#'
#' @returns For a single file: an object of class flowstate; for multiple files: a list of flowstate objects
#' @export
#'
#' @examples
#' fcs.file.paths <- system.file("extdata", package = "flowstate") |>
#' list.files(full.names = TRUE, pattern = "BLOCK.*.fcs")
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
#' #.fcs TEXT segment parsed and stored as three elements (data.tables):
#'   fs$parameters #instrument-specific measurement parameters
#'   fs$keywords #instrument/sample-specific metadata
#'   fs$spill #instrument/sample-specific spillover
#'
#' #read all .fcs files as flowstate objects; concatenate into a single object
#' fs <- read.flowstate(
#'   fcs.file.paths,
#'   colnames.type="S",
#'   cofactor = 5000,
#'   concatenate = TRUE
#' )
#' class(fs)
#' fs$keywords
#'
read.flowstate<-function(
    fcs.file.paths,
    colnames.type=c("S_N","S","N"),
    S.func=NULL,
    cofactor=5000,
    sample.id='TUBENAME',
    concatenate=FALSE
)
{
  ##add names to fcs.files.paths
  fcs.file.paths<-stats::setNames(fcs.file.paths,nm=basename(fcs.file.paths))
  ##create the object(s)
  fs<-lapply(fcs.file.paths,function(fcs.file.path){
    message(paste(basename(fcs.file.path),"-->","flowstate"))
    flowstate.from.file.path(fcs.file.path,S.func = S.func)
  })
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
    })#define old here as fs.obj$parameters[['N']]?
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
        ##Cytek Aurora; spectral keyword/value; $PnTYPE/Unmixed_Fluorescence
        cols.transform<-fs.obj$parameters[grep("[RawUnmixed]_Fluorescence",TYPE)][[colnames.type]]
        ##
        message(paste(fs.obj$keywords[,`$FIL`],"-->","transforming..."))
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
    if(concatenate){
      fs <- concatenate.flowstate(fs)
    }else{
      return(fs)
    }
  }
}
