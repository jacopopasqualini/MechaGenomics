import os

configfile: "../configuration/ecoGutConfig.json"

PARAMS=config["rules_parameters"]
THREADS=config["threads"]
DATABASES=config["databases"]

rule kaijuUHGG:
	input:
		fastqD="../data/decontaminatedFastq/{sample}_decon.fastq.gz"
		
	output:
		kaiju_out="../data/kaijuUHGG/{sample}.out.gz",
		kaiju_tab="../data/kaijuUHGG/{sample}_uhgg.tsv"

	params:
		uhgg=DATABASES["kaijuUHGG"]
		
	conda: 
		"../envs/metagenomics.yml"
	
	threads: 
		config["threads"]["kaijuUHGG"]
	
	message: 
		">>> Kaiju UHGG v2 Classification."

	log: 
		"../data/kaijuUHGG/log/{sample}-kaiju-uhgg.log"
	
	shell:
		""" kaiju -v -z 16 -t {params.uhgg}/nodes.dmp -f {params.uhgg}/UHGG.fmi -i {input.fastqD} -o ../data/kaijuUHGG/uc_{wildcards.sample}.out && """
		""" kaiju2table -t {params.uhgg}/nodes.dmp -n {params.uhgg}/names.dmp -r species -o {output.kaiju_tab} ../data/kaijuUHGG/uc_{wildcards.sample}.out && """
		""" grep -v 'U' ../data/kaijuUHGG/uc_{wildcards.sample}.out > ../data/kaijuUHGG/{wildcards.sample}.out && """
		""" gzip ../data/kaijuUHGG/{wildcards.sample}.out && """
		""" rm ../data/kaijuUHGG/uc_{wildcards.sample}.out """

rule kaijuRefSeq:
	input:
		fastqD="../data/decontaminatedFastq/{sample}_decon.fastq.gz"
		
	output:
		kaiju_out="../data/kaijuRefSeq/{sample}.out.gz",
		kaiju_tab="../data/kaijuRefSeq/{sample}_refseq.tsv"

	params:
		refseq=DATABASES["kaijuRefSeq"]
		
	conda: "../envs/metagenomics.yml"
	
	threads: 
		config["threads"]["kaijuRefSeq"]
	
	message:
		">>> Kaiju RefSeq Classification."

	log: 
		"../data/kaijuRefSeq/log/{sample}-kaiju-refseq.log"
	
	shell:
		""" kaiju -v -z 16 -t {params.refseq}/nodes.dmp -f {params.refseq}/kaiju_db_refseq.fmi -i {input.fastqD} -o ../data/kaijuRefSeq/uc_{wildcards.sample}.out && """
		""" kaiju2table -t {params.refseq}/nodes.dmp -n {params.refseq}/names.dmp -r species -o {output.kaiju_tab} ../data/kaijuRefSeq/uc_{wildcards.sample}.out && """
		""" grep -v 'U' ../data/kaijuRefSeq/uc_{wildcards.sample}.out > ../data/kaijuRefSeq/{wildcards.sample}.out && """
		""" gzip ../data/kaijuRefSeq/{wildcards.sample}.out && """
		""" rm ../data/kaijuRefSeq/uc_{wildcards.sample}.out """

rule kaijuCorePFAM:
	input:
		fastqD="../data/decontaminatedFastq/{sample}_decon.fastq.gz"
		
	output:
		kaiju_out="../data/kaijuCorePFAM/{sample}.out.gz",
		kaiju_tab="../data/kaijuCorePFAM/{sample}.tsv"

	params:
		core_pfam=DATABASES["kaijuCorePFAM"]
		
	conda: "../envs/metagenomics.yml"
	
	threads: 
		config["threads"]["kaijuCorePFAM"]
	
	message: 
		">>> Kaiju Core-PFAM Classification."

	log: "../data/kaijuCorePFAM/log/{sample}-kaiju-core.log"
	
	shell:
		""" kaiju -v -z 16 -t {params.core_pfam}/nodes.dmp -f {params.core_pfam}/corePFAM.fmi -i {input.fastqD} -o ../data/kaijuCorePFAM/uc_{wildcards.sample}.out && """
		""" kaiju2table -t {params.core_pfam}/nodes.dmp -n {params.core_pfam}/names.dmp -r species -o {output.kaiju_tab} ../data/kaijuCorePFAM/uc_{wildcards.sample}.out && """
		""" grep -v 'U' ../data/kaijuCorePFAM/uc_{wildcards.sample}.out > ../data/kaijuCorePFAM/{wildcards.sample}.out && """
		""" gzip ../data/kaijuCorePFAM/{wildcards.sample}.out && """
		""" rm ../data/kaijuCorePFAM/uc_{wildcards.sample}.out """

rule kaijuCoreUHGG:
	input:
		fastqD="../data/decontaminatedFastq/{sample}_decon.fastq.gz"
		
	output:
		kaiju_out="../data/kaijuCoreUHGG/{sample}.out.gz",
		kaiju_tab="../data/kaijuCoreUHGG/{sample}.tsv"

	params:
		core_uhgg=DATABASES["kaijuCoreUHGG"]
		
	conda: "../envs/metagenomics.yml"
	
	threads: 
		config["threads"]["kaijuCoreUHGG"]
	
	message: 
		">>> Kaiju Core-PFAM Classification based on UHGG MAGs."

	log: "../data/kaijuCore/log/{sample}-kaiju-core-uhgg.log"
	
	shell:
		""" kaiju -v -a mem -z 16 -t {params.core_uhgg}/nodes.dmp -f {params.core_uhgg}/coreUHGGv2.fmi -i {input.fastqD} -o ../data/kaijuCoreUHGG/uc_{wildcards.sample}.out && """
		""" kaiju2table -t {params.core_uhgg}/nodes.dmp -n {params.core_uhgg}/names.dmp -r species -o {output.kaiju_tab} ../data/kaijuCoreUHGG/uc_{wildcards.sample}.out && """
		""" grep -v 'U' ../data/kaijuCoreUHGG/uc_{wildcards.sample}.out > ../data/kaijuCoreUHGG/{wildcards.sample}.out && """
		""" gzip ../data/kaijuCoreUHGG/{wildcards.sample}.out && """
		""" rm ../data/kaijuCoreUHGG/uc_{wildcards.sample}.out """