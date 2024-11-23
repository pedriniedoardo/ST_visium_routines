library(GEOquery)
gse <- getGEO('GSE222322',GSEMatrix = F,destdir = "test/")

folder <- names(GSMList(gse))

lapply(folder, function(x){
  dir.create(paste0("test/",x))
})

gse@gsms$GSM6919905@header$supplementary_file_2

download.file(gse@gsms$GSM6919905@header$supplementary_file_1,destfile = paste0("test/GSM6919905/filtered_feature_bc_matrix.h5"))

pmap(list(gse@gsms,names(gse@gsms)),function(x,name){
  download.file(x@header$supplementary_file_1,destfile = paste0("test/",name,"/"))
  gse@gsms$GSM6919905@header$supplementary_file_1
})

file_name <- basename(gse@gsms$GSM6919905@header$supplementary_file_1) %>% 
  str_remove_all("_filtered_feature_bc_matrix.h5")

write_tsv(x = data.frame("test"),file = paste0("test/",name,"/",file_name,".txt"))

