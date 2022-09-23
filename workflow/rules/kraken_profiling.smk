configfile: "../configuration/ecoGutConfig.json"

import os

DATA_DIR=config["data_directories"]
SAMPLES=config["samples"]
PARAMS=config["rules_parameters"]
THREADS=config["threads"]
DATABASES=config["databases"]

rule krakenUHGG:
	input:
		community_fastq=expand("{DECONTAM_DIR}/{sample}_decon.fastq.gz",DECONTAM_DIR=DATA_DIR["decontaminatedFastq"],sample=SAMPLES),
		
	output:
		kraken=expand("{KR_UHGG_OUT}/{sample}.out.gz",KR_UHGG_OUT=DATA_DIR["krakenUHGGOut"],sample=SAMPLES),
		kraken_report=expand("{KR_UHGG_OUT}/{sample}_report.txt",KR_UHGG_OUT=DATA_DIR["krakenUHGGOut"],sample=SAMPLES)
	params:
		sample=expand("{sample}",sample=SAMPLES),
		OUT_DIR=DATA_DIR["krakenUHGGOut"],
		uhgg=DATABASES["krakenUHGG"],
		
	conda: "../envs/metagenomics.yml"
	
	threads: THREADS["krakenUHGG"]
	
	message: ">>> Kraken UHGG Classification."

	log: expand('{LOG_DIR}/{sample}-kraken-uhgg.log', LOG_DIR=os.path.join(DATA_DIR["krakenUHGGOut"],'log'), sample=SAMPLES)
	
	shell:
		""" kraken2 --gzip-compressed --db {params.uhgg} --output {params.OUT_DIR}/{params.sample}.out --report {output.kraken_report} {input.community_fastq} && """
		""" gzip {params.OUT_DIR}/{params.sample}.out """