# Prepare Analysis Project/Directories

## Create a data analysis project

``` r
##create an analysis project with a meaningful name -- here the name of the study
path.project <- file.path(tempdir(),"COVAIL")
usethis::create_project(
  path = path.project,
  open = FALSE
) |> suppressMessages()
##within the project: create a directory structure for organizing analysis/workflow
sapply(paste0("data_",c('source','modified','results','meta','temp')),function(.dir){
    .dir <- file.path(path.project,.dir)
    if(!dir.exists(.dir)) dir.create(.dir)
}) |> invisible()
##directory structure of project
list.dirs(path.project)
#> [1] "/tmp/Rtmpju6t7M/COVAIL"              
#> [2] "/tmp/Rtmpju6t7M/COVAIL/data_meta"    
#> [3] "/tmp/Rtmpju6t7M/COVAIL/data_modified"
#> [4] "/tmp/Rtmpju6t7M/COVAIL/data_results" 
#> [5] "/tmp/Rtmpju6t7M/COVAIL/data_source"  
#> [6] "/tmp/Rtmpju6t7M/COVAIL/data_temp"    
#> [7] "/tmp/Rtmpju6t7M/COVAIL/R"
```

## Source files

``` r
temp <- tempfile(tmpdir = file.path(path.project,"data_temp"),fileext = ".zip")
zip.link <- "https://figshare.com/ndownloader/articles/30380776?private_link=38afd2ef2906779e8605"
curl::curl_download(url = zip.link, destfile = temp)
##files contained in the .zip
utils::unzip(temp,list = TRUE) |> head()
##unzip all to "./data_source"
utils::unzip(temp,files=NULL,exdir = file.path(path.project,"data_source"))
##source files (raw) in "./data_source"
files.source <- list.files(file.path(path.project,"data_source"),full.names=T,recursive = T)
print(sub("^.*/COVAIL","",files.source)) |> head()
##verify integrity using hashes
hash.historic <- "315cb3f534f175fa070df488a49f525a"
hash.source <- sapply(files.source,tools::md5sum) |> as.vector() |> digest::digest()
hash.historic == hash.source
if(hash.historic == hash.source){
  unlink(temp)
}
```
