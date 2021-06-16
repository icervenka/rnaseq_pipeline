#!/usr/bin/Rscript

library(rmarkdown)
library(rjson)

#set working directory
setwd(directory)

# import script variables from json -------------------------------------------------------------
global_vars <- fromJSON(file="global_vars.json")

# set names of files to be copied and used
from = "/usr/local/bin/r_rnaseq"
from_rmd_child = paste0(from,"/rmd_child")
to = directory
to_rmd_child = paste0(to,"/rmd_child")

report_name = paste0("report_", global_vars$file_prefix, ".Rmd")
rmd_files = list.files(from_rmd_child, full.names = TRUE)
rmd_main = paste0(from,"/",report_name)

# copy files from script location to local directory
dir.create(to_rmd_child)
file.copy(rmd_files, to_rmd_child)
file.copy(rmd_main, to)

# render html document
render(report_name, 
       output_dir = paste0(global_vars$output_dir), 
       output_format = "html_document")

# remove files from script directory
rmd_del_files = list.files(to_rmd_child, full.names = TRUE)
file.remove(rmd_del_files)
rmd_del_main = paste0(to,"/",report_name)
file.remove(rmd_del_main)
res = unlink(to_rmd_child, force = TRUE)
