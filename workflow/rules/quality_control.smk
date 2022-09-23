import os

configfile: "../configuration/ecoGutConfig.json"


rule trimmedFastqQC:
	input:
		fastqD="../data/decontaminatedFastq/{sample}_decon.fastq.gz"

	output:
		htmlD="../data/decontaminatedFastQC/{sample}_decon_fastqc.html",
		zipD="../data/decontaminatedFastQC/{sample}_decon_fastqc.zip"

	conda: 
		"../envs/metagenomics.yml"

	threads: 
		config["threads"]["rawFastqQC"]

	log:
		"../data/decontaminatedFastQC/log/{sample}-qc-trim.log"

	message: 
		">>> Quality assesment of raw fastq files AFTER trimming and decontamination."	

	shell:
		"fastqc  --outdir ../data/decontaminatedFastQC {input.fastqD}"

rule decontaminated2table:
	input:
		fastqD="../data/decontaminatedFastq/{sample}_decon.fastq.gz"

	output:
		deconTab="../data/decontaminatedTable/{sample}_decon_fx2tab.out"

	conda: 
		"../envs/metagenomics.yml"

	threads:
		config["threads"]["fastq2table"]

	log:
		"..data/decontaminatedTable/log/{sample}-fx2tab.log"

	message:
		">>> Quality resume of raw fastq files AFTER trimming and decontamination."

	shell:
		"seqkit fx2tab {input.fastqD} --name --avg-qual --length --out-file {output.deconTab}"