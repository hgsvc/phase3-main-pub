rule C01_linkGeneModel:
    input:
        gm=config["genemodel"]
    output:
        gmlink="results/C01_gene_model.fasta"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/C01_linkGeneModel.log"
    shell:"""
        echo "##### C01_linkGeneModel" > {log}
        ln -s {input.gm} {params.workflowDir}/../{output.gmlink} 2>> {log}
    """

rule C02_SimplifyFastaGeneDefs:
    input:
        gm="results/C01_gene_model.fasta"
    output:
        gene_names_key_noSuffix=temp("results/C02_gene_names_key_noSuffix.tsv"),
        gene_names_key="results/C02_gene_names_key.tsv",
        gm_simp="results/C02_gene_model_renamed.fasta"
    params:
        assumeClearAndUniqueCodes=config["flag_assume_clear_and_unique_gene_codes"],
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/C02_SimplifyFastaGeneDefs.log"
    shell:"""
    {{
        # Assumes Gene Code is in last set of parentheses in fasta entry descriptor.
        echo "##### C02_SimplifyFastaGeneDefs" > {log}
        if [ {params.assumeClearAndUniqueCodes} = "True" ]
        then
            echo "### Do not Simplify Gene Model Gene Defs: create symlink instead." >> {log}
            ln -s {params.workflowDir}/../{input.gm} {params.workflowDir}/../{output.gm_simp}
            touch {output.gene_names_key_noSuffix}
            touch {output.gene_names_key}
        else
            echo "### Parse Gene Code" >> {log}
            cat {input.gm} | awk 'BEGIN {{OFS="\\t"}} \
                (/^>/)\
                    {{split($0,openPieces,"("); \
                    split(openPieces[length(openPieces)],closedPiece,")"); \
                    geneCode=closedPiece[1]; \
                    print NR,geneCode,$0}}' 1> {output.gene_names_key_noSuffix}

            echo "### Add Suffix to Gene Code" >> {log}
            cat {output.gene_names_key_noSuffix} | \
                sort -k2,2 -k1,1n | \
                awk 'BEGIN {{FS="\\t"; OFS="\\t"; lastCode=""; geneCount=1}} \
                    {{if ($2==lastCode) \
                        {{geneCount++}} \
                    else \
                        {{geneCount=1}} \
                    print $1,$2"|"geneCount,$3; \
                    lastCode=$2}}' | \
                sort -k1,1n 1> {output.gene_names_key}

            echo "### Replace fasta descriptors with Gene Code" >> {log}
            cat {input.gm} | \
                awk ' \
                    !(/^>/) \
                        {{print $0}} \
                    (/^>/) \
                        {{getline < "{output.gene_names_key}"; \
                        print ">"$2}}' 1> {output.gm_simp}
        fi
    }} 2>> {log}
    """

rule C03_RemoveUncharacterizedGenes:
    input:
        gm="results/C02_gene_model_renamed.fasta"
    output:
        gm_characterized="results/C03_gene_model_filt.fasta"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/C03_RemoveUncharacterizedGenes.log"
    params:
        filt_uncharacterized=config["flag_filt_uncharacterized_genes"],
        uncharacterized_prefix=config["uncharacterized_gene_name_prefix"],
        workflowDir=workflow.basedir
    shell:"""
    {{
        echo "##### C03_RemoveUncharacterizedGenes" > {log}
        if [ {params.filt_uncharacterized} = "True" ]
        then
            echo "### Remove Uncharacterized Genes" >> {log}
            cat {input.gm} | awk -v unchar_prefix={params.uncharacterized_prefix} ' \
                (/^>/) \
                    {{ prefix=substr($0,1,length(unchar_prefix)+1) }} \
                (prefix!=">"unchar_prefix) \
                    {{ print $0 }}' 1> {output.gm_characterized}
        else
            echo "### Do not remove uncharacterized genes: create symlink instead." >> {log}
            ln -s {params.workflowDir}/../{input.gm} {params.workflowDir}/../{output.gm_characterized}
        fi
    }} 2>> {log}
    """
