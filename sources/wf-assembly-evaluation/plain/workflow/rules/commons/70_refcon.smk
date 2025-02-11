if USE_REFERENCE_CONTAINER:

    localrules:
        refcon_dump_manifest,
        refcon_cache_manifests,

    rule refcon_dump_manifest:
        input:
            sif=DIR_REFCON.joinpath("{refcon_name}.sif"),
        output:
            manifest=DIR_PROC.joinpath(".cache", "refcon", "{refcon_name}.manifest"),
        envmodules:
            ENV_MODULE_SINGULARITY,
        shell:
            "{input.sif} manifest > {output.manifest}"

    rule refcon_run_get_file:
        """
        Snakemake interacts with Singularity containers using "exec",
        which leads to a problem for the "refcon_run_get_file".
        Dynamically setting the Singularity container for the
        "container:" (former "singularity") keyword results in a parsing error for
        unclear reasons. Hence, for now, force the use of
        "singularity run" to extract data from reference containers
        (i.e., treat them like a regular file)
        """
        input:
            cache=trigger_refcon_manifest_caching,
        output:
            DIR_GLOBAL_REF.joinpath("{filename}"),
        envmodules:
            ENV_MODULE_SINGULARITY,
        params:
            refcon_path=lambda wildcards, input: refcon_find_container(
                input.cache, wildcards.filename
            ),
            acc_ref=lambda wildcards, output: register_reference(output),
        shell:
            "{params.refcon_path} get {wildcards.filename} {output}"

    checkpoint refcon_cache_manifests:
        input:
            manifests=expand(
                DIR_PROC.joinpath(".cache", "refcon", "{refcon_name}.manifest"),
                refcon_name=load_reference_container_names(),
            ),
        output:
            cache=DIR_PROC.joinpath(".cache", "refcon", "refcon_manifests.cache"),
        run:
            import pandas

            merged_manifests = []
            for manifest_file in input.manifests:
                container_name = pathlib.Path(manifest_file).stem
                manifest = pandas.read_csv(manifest_file, sep="\t", header=0)
                manifest["refcon_name"] = container_name
                merged_manifests.append(manifest)
            merged_manifests = pandas.concat(
                merged_manifests, axis=0, ignore_index=False
            )

            merged_manifests.to_csv(
                output.cache, header=True, index=False, sep="\t"
            )
            # END OF RUN BLOCK
