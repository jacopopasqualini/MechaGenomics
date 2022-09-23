configfile: "../configuration/ecoGutConfig.json"

import os

rule hostDecontamination:
	input:
		fastqT="../data/trimmedFastq/{sample}_trimmed.fastq.gz"
		
	output:
		fastqK="../data/decontaminatedFastq/{sample}_trimmed_kneaddata.fastq",
		fastqC="../data/decontaminatedFastq/{sample}_trimmed_kneaddata_hg37dec_v0.1_bowtie2_contam.fastq",
		deconLog="../data/decontaminatedFastq/log/{sample}_decontamination.log"

	params:
		hg37DB=config["databases"]["hg37BowtieDB"],
		bowtie2_mode = config["rules_parameters"]["hostDecontamination"]["bowtie2_mode"]
		
	conda: 
		"../envs/metagenomics.yml"
	
	threads: 
		config["threads"]["hostDecontamination"]
	
	message: 
		">>> Human host decontamination."

	log: 
		"../data/decontaminatedFastq/log/{sample}-decontaminationKneaddata.log"
	
	shell:
		""" kneaddata --input {input.fastqT} \
		--output ../data/decontaminatedFastq --reference-db {params.hg37DB} \
		--bypass-trf --bypass-trim --remove-intermediate-output \
		--threads {threads} --bowtie2-options="{params.bowtie2_mode}" \
		--log {output.deconLog}
		"""

rule cleanDecontamination:
	input:
		fastqK="../data/decontaminatedFastq/{sample}_trimmed_kneaddata.fastq",
		fastqC="../data/decontaminatedFastq/{sample}_trimmed_kneaddata_hg37dec_v0.1_bowtie2_contam.fastq",
		fastqT="../data/trimmedFastq/{sample}_trimmed.fastq.gz"

	output:
		fastqD="../data/decontaminatedFastq/{sample}_decon.fastq.gz"
	
	message:
		">>> Cleaning decontamination session."

	log:
		"../data/decontaminatedFastq/log/{sample}-clean-decontamination.log"
	
	run:
		shell("gzip -c {input.fastqK} > {output.fastqD}")
		shell("rm {input.fastqK} {input.fastqC} {input.fastqT}")