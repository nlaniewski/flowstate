# dir.fcs<-("W:/COVAIL_Rochester23X075F/data_source/COVAIL_002_CYTOKINE_2025-02-27/Unmixed/BarcodedPool/")
# fcs.file.paths<-list.files(dir.fcs,full.names = T,pattern = ".fcs")
#
# invisible(
#   lapply(grep("BLOCK1",fcs.file.paths,value = T),function(fcs.file.path){
#     ##read .fcs file; no transformation
#     fs<-flowstate::read.flowstate(fcs.file.path,cofactor = NULL)
#     ##subset [['data']]; first 2E3 rows
#     fs$data<-fs$data[1:2E3]
#     ##write out the subset
#     flowstate:::write.flowstate(
#       fs,
#       filename = file.path("W:/fcs_test_files/",fs$keywords$`$FIL`)
#     )
#   })
# )
