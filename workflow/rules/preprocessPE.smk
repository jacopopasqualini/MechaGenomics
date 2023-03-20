configfile: "../configuration/ecoGutConfig.json"

import os
import glob

PARAMS=config["rules_parameters"]
THREADS=config["threads"]

rule rawFastqQC:
	input:
		fastqF="../data/rawFastq/{sample}_1.fastq.gz",
		fastqR="../data/rawFastq/{sample}_2.fastq.gz"
	output:
		htmlF="../data/rawFastQC/{sample}_1_fastqc.html",
		zipF="../data/rawFastQC/{sample}_1_fastqc.zip",
		htmlR="../data/rawFastQC/{sample}_2_fastqc.html",
		zipR="../data/rawFastQC/{sample}_2_fastqc.zip"
	conda: 
		"../envs/metagenomics.yml"
	threads: 
		THREADS["rawFastqQC"]
	log: 
		"../data/rawFastQC/log/{sample}-qc-before-trim.log"
	message: 
		">>> Quality assesment of raw fastq files BEFORE trimming and decontamination."
	shell:
		"""fastqc  --outdir ../data/rawFastQC {input.fastqF} && """
		"""fastqc  --outdir ../data/rawFastQC {input.fastqR} """


rule fastqTrimming:
	input:
		fastq1="../data/rawFastq/{sample}_1.fastq.gz",
		fastq2="../data/rawFastq/{sample}_2.fastq.gz"
		
	output:
		fastqM="../data/trimmedFastq/{sample}_merged.fastq.gz",
		fastq1P="../data/trimmedFastq/{sample}_1P.fastq.gz",
		fastq2P="../data/trimmedFastq/{sample}_2P.fastq.gz",
		fastq1U="../data/trimmedFastq/{sample}_1U.fastq.gz",
		fastq2U="../data/trimmedFastq/{sample}_2U.fastq.gz",
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
		""" fastp --in1 {input.fastq1} --in2 {input.fastq2} \
			--length_required {params.reads_length} --average_qual {params.avg_phred} \
			--merge --merged_out {output.fastqM} \
			--out1  {output.fastq1P} --out2  {output.fastq2P} \
			--unpaired1  {output.fastq1U} --unpaired2  {output.fastq2U} \
			--json {output.fastpJson} \
			--html {output.fastpHtml} """


rule cleanTrimming:
	input:
		fastqM="../data/trimmedFastq/{sample}_merged.fastq.gz",
		fastq1P="../data/trimmedFastq/{sample}_1P.fastq.gz",
		fastq2P="../data/trimmedFastq/{sample}_2P.fastq.gz",
		fastq1U="../data/trimmedFastq/{sample}_1U.fastq.gz",
		fastq2U="../data/trimmedFastq/{sample}_2U.fastq.gz"

	output:
		fastqT="../data/trimmedFastq/{sample}_trimmed.fastq.gz"

	conda: 
		"../envs/metagenomics.yml"

	message:
		">>> Clean trimming session."
	log: 
		"../data/trimmedFastq/log/{sample}-clean-trim.log"

	shell:
		"seqkit concat {input.fastq1P} {input.fastq2P} > ../data/trimmedFastq/{wildcards.sample}_concat.fastq & "
		"gzip -c ../data/trimmedFastq/{wildcards.sample}_concat.fastq >  ../data/trimmedFastq/{wildcards.sample}_concat.fastq.gz & "
		"cat {input.fastqM} ../data/trimmedFastq/{wildcards.sample}_concat.fastq.gz > {output.fastqT} & "
		"rm {input.fastq1U} {input.fastq2U} {input.fastq1P} {input.fastq2P} {input.fastqM} & "
		"rm ../data/trimmedFastq/{wildcards.sample}_concat.fastq ../data/trimmedFastq/{wildcards.sample}_concat.fastq.gz"
