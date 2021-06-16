#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 16:49:31 2017

@author: icervenka
"""

core_no = 16

genome_index_path_m = "/shared/genome/GRCm38/hisat/GRCm38"
annotation_path_m = "/shared/genome/GRCm38/annotation/GRCm38.gtf"
annotation_path_flat_m = "/shared/genome/GRCm38/annotation/GRCm38_flat.gff"

genome_index_path_h = "/shared/genome/GRCh38/hisat/GRCh38"
annotation_path_h = "/shared/genome/GRCh38/annotation/GRCh38.gtf"
annotation_path_flat_h = "/shared/genome/GRCh38/annotation/GRCh38_flat.gff"

hisat_log_dir = ""
hisat_log_filename = ""
htseq_log_dir = ""
htseq_log_filename = ""

accepted_extensions = ["fa", "fq", "fastq", "fasta"]

htseq_pheno_basename = "htseq_pheno.txt"
R_path = "/usr/lib/R/bin/Rscript"
exec_path = "/usr/local/bin/r_rnaseq"
python_path = ""

