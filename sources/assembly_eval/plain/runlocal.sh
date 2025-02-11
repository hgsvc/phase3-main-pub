#!/usr/bin/env bash
set -euo pipefail


SOURCE="${BASH_SOURCE[0]}"

while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink                                                              
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink
done

DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


module load miniconda/23.5.2

snakefile=$DIR/assembly_eval.smk
waitTime=10
retry=1


logDir=log
mkdir -p $logDir

#
# run snakemake
#

snakemake -p \
        -s $snakefile \
        --latency-wait $waitTime \
        --use-singularity \
        --singularity-args="--bind /net/:/net/" \
        --restart-times $retry  \
        --ri \
        --jobs $@
