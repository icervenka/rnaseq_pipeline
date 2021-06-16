#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 11:49:59 2017

@author: icervenka
"""

from collections import Counter
from subprocess import call
import pandas as pd
import glob
import os
import os.path


def _itersplit(l, splitters):
    current = []
    for item in l:
        if item in splitters:
            yield current
            current = []
        else:
            current.append(item)
    yield current


def listsplit(l, *splitters):
       return [subl for subl in _itersplit(l, splitters) if subl] 
  
    
def compare(s, t):
    return(Counter(s) == Counter(t))


def base_file_vect(arr, add_cleanup = []):
    ret_arr = [os.path.basename(x) for x in arr]
    ret_arr = [os.path.splitext(x)[0] for x in ret_arr]
    for item in add_cleanup:
        ret_arr = [x.replace(item, "") for x in ret_arr]
    return(ret_arr)


def infer_extension(extensions, directory):
    num_files = {}
    for ext in extensions:
        file_count = len(glob.glob(directory+"/*."+ext))
        num_files[ext] = file_count
    
    non_null = [k for k, v in num_files.items() if v > 0]
    if len(non_null) != 1:
        print("Different file types present in folder, unable to choose.")
    else:
        return(non_null[0])


def parse_hisat_log(basedir = ".", filePath = "hisat_log.txt"):
    filelines = None
    filename_list = []
    alignment_list = []
    columns = ['file', 'total','align_zero','align_once','align_more']
    
    with open(basedir+"/"+filePath) as file:
        filelines = file.read().splitlines()
    
    filelines = [x.strip() for x in filelines]
    log_list = listsplit(filelines, '')
    
    for item in log_list:
        filename_list.append(item[0])
        al = []
        al = [item[1], item[3], item[4], item[5]]
        al = [x.split(" ", 1)[0] for x in al]  
        alignment_list.append(al)
    
    filename_list = base_file_vect(filename_list)
    df_hisatlog = pd.concat([pd.Series(filename_list), pd.DataFrame(alignment_list)], axis=1)
    df_hisatlog.columns = columns
    
    #TODO add percentages of aligned reads to dataframe
    return(df_hisatlog)


def parse_htseq_log(htseq_files, basedir = "."):
    htseq_log_arr = []
    column_names = ['file', 'no_feature', 'ambiguous', 'low_quality', 'not_aligned', 'not_unique']
    #files = glob.glob(basedir+"/"+"*.txt")
    
    for file in htseq_files:
        with open(file+".txt") as f: 
            file_lines = f.read().splitlines()
            log_array = file_lines[-5:] 
            vals = []
            for item in log_array:
                val = item.split("\t")[1]
                vals.append(val)       
        htseq_log_arr.append(vals)   
            
    df_htseq = pd.DataFrame(htseq_log_arr)
    files = htseq_files
    
    df_htseq = pd.concat([pd.Series(files), df_htseq], axis=1)
    df_htseq.columns = column_names
    df_htseq = df_htseq.set_index('file')
    df_htseq = df_htseq.astype(int)
    df_htseq = df_htseq.reset_index()
    
    return(df_htseq)


def pheno_exists(directory, filename = 'htseq_pheno.txt'):
    os.chdir(directory)
    if not os.path.isfile(filename):
        filename = input("Please provide name of existing sample mapping file. File {0} does not exists.".format(filename))
        return(filename)
    else:
        return(filename)
    

def verify_pheno(num_samples, ref_level, num_contrasts = 1, directory = ".", filename = "htseq_pheno.txt"):
        
    htseq_pheno = pd.read_csv("htseq_pheno.txt", sep = ",", header = None)
    
    conditions = htseq_pheno.iloc[:,2]
    
    if not conditions.str.contains(ref_level).any():
        print("Specified reference level is not present in conditions")
        return(False)
        
    num_rows = htseq_pheno.shape[0]
    num_columns = htseq_pheno.shape[1]
    
    if num_samples < num_rows:
        print("Number of rows in mapping file doesn't correspond to number of "+
        "samples. Please provide {0} lines in mapping file".format(num_samples))
        return(False)
    
    if num_columns != (2 + num_contrasts):
        print("Number of columns does not match. Please provide column for:\n" +
        "1) sample name, \n2) file name \n3) {0} contrasts".format(num_contrasts))
        return(False)
    
    return(True)
    

        
def verify_deg(directory, run = 'h'):
    files = [x for x in glob.glob(directory+"/*.txt")]
    if files:
        return(True)
    else:
        print("No files with counts of genome features to process")
        return(False)
    
    
def verify_markdown(directory, r_script_prefix):     
    if not os.path.isfile(directory+"/"+r_script_prefix+".RData"):
        print("No correct RData file is present.")
    else:
        return(True)
        
        
def run_command(x):
    call(x.split(), shell=False)