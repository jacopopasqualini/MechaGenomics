import os

configfile: "../configuration/ecoGutConfig.json"

PARAMS=config["rules_parameters"]
THREADS=config["threads"]
DATABASES=config["databases"]

rule krakenUHGG:
	input:
		community_fastq="../data/decontaminatedFastq/{sample}_decon.fastq.gz"
		
	output:
		kraken_out="../data/krakenUHGG/{sample}.out.gz",
		kraken_report="../data/krakenUHGG/{sample}.report"
	params:
		uhgg=DATABASES["krakenUHGG"],
		
	conda: "../envs/metagenomics.yml"
	
	threads: THREADS["krakenUHGG"]
	
	message: ">>> Kraken UHGG Classification."

	log: 
		"../data/krakenUHGG/log/{sample}-KrakenUHGG.log"
	
	shell:
		""" kraken2 --gzip-compressed --db {params.uhgg} --output ../data/krakenUHGG/{wildcards.sample}.out --report {output.kraken_report} {input.community_fastq} && """
		""" gzip -c ../data/krakenUHGG/{wildcards.sample}.out > {output.kraken_out} &&
		    rm ../data/krakenUHGG/{wildcards.sample}.out """

rule krakenCoreUHGG:
	input:
		community_fastq="../data/decontaminatedFastq/{sample}_decon.fastq.gz"
		
	output:
		kraken_out="../data/krakenCoreUHGG/{sample}.out.gz",
		kraken_report="../data/krakenCoreUHGG/{sample}.report"
	params:
		uhgg=DATABASES["krakenCoreUHGG"],
		
	conda: "../envs/metagenomics.yml"
	
	threads: THREADS["krakenCoreUHGG"]
	
	message: ">>> Kraken CoreUHGG Classification."

	log: 
		"../data/krakenUHGG/log/{sample}-KrakenCoreUHGG.log"
	
	shell:
		""" kraken2 --gzip-compressed --db {params.uhgg} --output ../data/krakenCoreUHGG/{wildcards.sample}.out --report {output.kraken_report} {input.community_fastq} && """
		""" gzip -c ../data/krakenCoreUHGG/{wildcards.sample}.out > {output.kraken_out} &&
		    rm ../data/krakenCoreUHGG/{wildcards.sample}.out """