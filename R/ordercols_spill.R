# fs<-read.flowstate(
#   "W:/COVAIL_Rochester23X075F/data_source/COVAIL_001_CYTOKINE_2025-02-13/Unmixed/BarcodedPool/COVAIL_001_CYTOKINE_BCpool_BLOCK1_1.fcs",
#   colnames.type = "S"
# )[[1]]
#
# data.table::setorder(
#   x = fs$parameters[
#     ,
#     ord := match(S.alias,stringr::str_sort(S.alias,numeric = T))
#   ],
#   "ord"
# )[,ord := NULL]
# fs$parameters[,par := paste0("$P",seq(.N))]
#
# data.table::setcolorder(
#   x = fs$data,
#   neworder = fs$parameters[['S.alias']]
# )
#
# set.seed(1337)
# plot.flowstate(fs,CD4,CD8,sample.n=5E4)
#
# str(fs$spill)
# spillover.update.value(fs,CD8,CD4,0.03)
# spillover.apply(fs)
# str(fs$spill)
#
# set.seed(1337)
# plot.flowstate(fs,CD4,CD8,sample.n=5E4)
#
# spillover.apply(fs,invert = TRUE)
# str(fs$spill)
#
# set.seed(1337)
# plot.flowstate(fs,CD4,CD8,sample.n=5E4)
