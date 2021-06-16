#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 11:49:59 2017

@author: icervenka
"""
#TODO remove current directory dots from helper function 
#TODO incorporate output dir in case current dir is not writable
#TODO hisat only run
#TODO flattened annotataion file check or generation
#TODO call R script with commandline arguments
#TODO hisat should accept variable log file name
#TODO change bash scripts into python scripts

from subprocess import call
import multiprocessing as mp
import os
import glob
import argparse
import json
import rnaseq_functions
import config

optParser = argparse.ArgumentParser(

   usage = "r_rnaseq [options] <directory> <ref_level>",

   description=
      "",

   epilog =
      "" )

optParser.add_argument("-p", "--paired", action="store_true",
   help = "Paired sequences")

optParser.add_argument( "-s", "--stranded", action="store_true",
   help = "Stranded protocol for RNA preparation")

optParser.add_argument('-l', "--align", action='store_true',
   help = "Align reads from sequencing files to target genome")

optParser.add_argument('-d', "--deg", action='store_true',
   help = "Perform analysis of differential gene expression")

optParser.add_argument( "-a", "--animal", type=str,
  choices = ( "m", "h", "n"), default = "m",
  help = "'m' - mouse or 'h' - human, 'n' - none, supply genome and annotation, default: mouse")

optParser.add_argument( "-c", "--count", type=str, dest="runtype",
  choices = ( "h", "s", "x", "n"), default = "h",
  help = "Type of analysis: 'h' - htseq-DESeq, 's' - stringtie-ballgown, 'x' - DEXSeq, 'n' - None, default: h")

optParser.add_argument( "-f", "--fdr", type=float, dest="fdr_value",
  default = 0.05,
  help = "FDR threshold value for significance of differentially expressed genes, default = 0.05")

optParser.add_argument("--index", type=str,
  help = "Genome index file for alignment. Required if animal is not specified")

optParser.add_argument("--annotation", type=str,
  help = "Annotation file for counting genomic features. Required if animal is not specified")

optParser.add_argument("directory", type=str, 
   help="Directory path to the source files")

optParser.add_argument("ref_level", type=str, 
   help="Name of reference group from htseq pheno file (control group name)")


# argument parsing and checking -----------------------------------------------
args = optParser.parse_args()

paired = 'yes' if args.paired else 'no'
stranded = 'yes' if args.stranded else 'no'
align = args.align
deg = args.deg
animal = args.animal
runtype = args.runtype
fdr_value = args.fdr_value
genome_index = args.index
annotation_file = args.annotation
directory = args.directory
ref_level = args.ref_level


if directory == ".":
    directory = os.getcwd()
elif directory.startswith('/'):
    pass
else:
    directory = os.getcwd()+"/"+directory

os.chdir(directory)


if animal == "m":
    genome_index = genome_index if genome_index else config.genome_index_path_m
    annotation_file = annotation_file if annotation_file else config.annotation_path_m
    flattened_annotation_file = config.annotation_path_flat_m
elif animal == "h":
    genome_index = genome_index if genome_index else config.genome_index_path_h
    annotation_file = annotation_file if annotation_file else config.annotation_path_h
    flattened_annotation_file = config.annotation_path_flat_h
else:
    if len(genome_index) < 1:
        print("Genome index was not specified")
    if len(annotation_file) < 1:
        print("Annotation file was not specified")
        

# files to process ------------------------------------------------------------
files = []

extension = rnaseq_functions.infer_extension(config.accepted_extensions, directory)
for file in glob.glob("*."+extension):
    files.append(file)

files_basename = [os.path.splitext(file)[0] for file in files]



# alignment to genome ---------------------------------------------------------
if align:    
    for file in files_basename:
        command = "r_hisat {0} {1}".format(directory+"/"+file, genome_index)
        call(command.split(), shell=False)

    df_hisatlog = rnaseq_functions.parse_hisat_log(basedir = directory)
    with open("hisat_log.csv", "w+") as file:
        df_hisatlog.to_csv(file, header = True, index = False)
 
       
# counting of genomic features ------------------------------------------------
if runtype == "h":
    r_script_prefix = "deseq"
    result_dir = directory+"/"+r_script_prefix
    
    htseq_pheno = rnaseq_functions.pheno_exists(directory)
    if not rnaseq_functions.verify_pheno(len(files_basename), ref_level):
        print("Problems with pheno file, please check.") # raise exception
    
    jobs = []
    for file in files_basename:
        command = "r_htseq {0} {1} {2}".format(directory+"/"+file, annotation_file, stranded)
        jobs.append(command)
    
    pool = mp.Pool(processes=(mp.cpu_count() - 1))
    pool.map_async(rnaseq_functions.run_command, jobs)
    pool.close()
    pool.join()

    df_htseqlog = rnaseq_functions.parse_htseq_log(files_basename, basedir = directory)   
    with open("htseq_log.csv", "w+") as file:
        df_htseqlog.to_csv(file, header = True, index = False)
       
elif runtype == "s":
    #call(["mkdir", "-p", "stringtie"])
    #call(["rm", "stringtie/mergelist.txt"])
    
    jobs = []
    for file in files_basename:
        command = "r_stringtie {0} {1} {2}".format(file, directory, annotation_file)
        jobs.append(command)
        #call(command.split(), shell=False)
     
    
    pool = mp.Pool(processes=(mp.cpu_count() - 1))
    pool.map_async(rnaseq_functions.run_command, jobs)
    pool.close()
    pool.join()
    
    call(["r_stringtiem", annotation_file], shell=False) 


elif runtype == "x":
    r_script_prefix = "dexeq"
    result_dir = directory+"/"+r_script_prefix
    
#    htseq_pheno = rnaseq_functions.pheno_exists(directory)
#    if not rnaseq_functions.verify_pheno(len(files_basename), ref_level):
#        print("Problems with pheno file, please check.") # raise exception
    
    jobs = []
    for file in files_basename:
        command = "r_dexseq {0} {1} {2} {3}".format(directory+"/"+file, flattened_annotation_file, stranded, paired)
        jobs.append(command)
    
    pool = mp.Pool(processes=(mp.cpu_count() - 1))
    pool.map_async(rnaseq_functions.run_command, jobs)
    pool.close()
    pool.join()
    
    df_dexseqlog = rnaseq_functions.parse_dexseq_log(files_basename, basedir = directory)   
    with open("dexseq_log.csv", "w+") as file:
        df_dexseqlog.to_csv(file, header = True, index = False)
     

else:
    r_script_prefix = "deseq"
    result_dir = directory+"/"+r_script_prefix
    htseq_pheno = rnaseq_functions.pheno_exists(directory)
    if not rnaseq_functions.verify_pheno(len(files_basename), ref_level):
        print("Problems with pheno file, please check.") # raise exception


    
# differential gene expression analysis ---------------------------------------
if deg and rnaseq_functions.verify_deg(directory):
    
    if not r_script_prefix:
        deg_package = ""
        while(deg_package not in ['deseq', 'stringtie', 'dexseq']):
            deg_package = input("Which DEG package do you wish to use (deseq/stringtie/dexseq):")
        r_script_prefix = deg_package
        result_dir = directory+"/"+r_script_prefix
            
    
    R_deg_args = {
    'directory':directory,
    'animal': animal,
    'fdr_value': fdr_value,
    'ref_level': ref_level,
    'file_prefix': r_script_prefix,
    'result_dir': result_dir,
    'pheno_file': htseq_pheno,
    'exec_path': config.exec_path
    }

    with open("global_vars.json", 'w') as jsonfile:
        json.dump(R_deg_args, jsonfile)
    
    rnaseq_args = [directory]
    call([config.R_path, config.exec_path+"/"+r_script_prefix+".R"] + rnaseq_args)



    