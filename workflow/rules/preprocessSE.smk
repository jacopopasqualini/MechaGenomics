configfile: "../configuration/ecoGutConfig.json"

import os

DATA_DIR=config["data_directories"]
SAMPLES=config["samples"]
PARAMS=config["rules_parameters"]
THREADS=config["threads"]

rule rawFastqQC:
	input:
		fastq="../data/rawFastq/{sample}_1.fastq.gz"

	output:
		html="../data/rawFastQC/{sample}_fastqc.html",
		zipf="../data/rawFastQC/{sample}_fastqc.zip"

	params: 
		OUT_DIR=DATA_DIR["rawFastQC"]

	conda: 
		"../envs/metagenomics.yml"

	threads: 
		THREADS["rawFastqQC"]

	log:
		"../data/rawFastQC/log/{sample}-qc-before-trim.log"
	
	message:
		">>> Quality assesment of raw fastq files BEFORE trimming and decontamination."
	
	shell:
		""" fastqc  --outdir ../data/rawFastQC/ {input.fastq} """

rule fastqTrimming:
	input:
		fastq="../data/rawFastq/{sample}.fastq.gz"
		
	output:
		fastqT="../data/trimmedFastq/{sample}_trimmed.fastq.gz",
		fastpJson="../data/fastpReport/{sample}_fastp.json",
		fastpHtml="../data/fastpReport/{sample}_fastp.html"

	params:
		reads_length=PARAMS["fastqTrimming"]["readsLength"],
		avg_phred=PARAMS["fastqTrimming"]["avgPhred"]

	conda: 
		"../envs/metagenomics.yml"
	
	threads: 
		THREADS["fastqTrimming"]
	
	message: 
		">>> Raw reads trimming."

	log: 
		"../data/trimmedFastq/log/{sample}-fastp.log"
	
	shell:
		""" fastp --in1 {input.fastq} \
		--length_required {params.reads_length} --average_qual {params.avg_phred} \
		--out1 {out.fastqT} --json {out.fastpJson} --html {out.fastpHtml} """

