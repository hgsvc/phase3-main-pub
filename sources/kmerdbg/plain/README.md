# Scripts to build and analyze compacted de Bruijn graphs

## Compacted de Bruijn graphs

We computed compacted de Bruijn graphs for all 3,202 samples from the [1000 Genomes high coverage project](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/). The pipeline is detailed below.

### CRAM to FASTQ conversion

We converted all CRAM files to FASTQ using [samtools](https://github.com/samtools/samtools).

```
samtools sort -n -o ${ID}.bam ${ID}.final.cram
samtools fastq -1 ${ID}.1.fq.gz -2 ${ID}.2.fq.gz -0 /dev/null -s /dev/null -n ${ID}.bam
```

### Kmer-based error correction

We next used [Lighter](https://github.com/mourisl/Lighter) for kmer-based error correction with a k-mer size of 23.

`lighter -r ${ID}.1.fq.gz -r ${ID}.2.fq.gz -trim -discard -k 23 3100000000 0.188`

### Compacted de Bruijn graph construction from error-corrected reads

[BCALM 2](https://github.com/GATB/bcalm) was used to compute the compacted de Bruijn graphs from the kmer-based error corrected reads.

```
ls -1 ${ID}.*.cor.fq.gz > ${ID}.list_reads
bcalm -in ${ID}.list_reads -kmer-size 61 -abundance-min 3
gzip ${ID}.unitigs.fa
```

The output of BCALM 2 is a set of unitigs (a non-branching path) of the de Bruijn graph. The unitig output format is described [here](https://github.com/GATB/bcalm#output). All 3,202 unitigs have been deposited at the [IGSR](https://www.internationalgenome.org/) FTP site in the [dbg_1kgp directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.0/dbg_1kgp/).

## Analysis of non-reference k-mers

We used [dicey](https://github.com/gear-genomics/dicey) and the [hs37d5](https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz) reference (GRCh37 coordinates) of the 1000 Genomes project to extract non-reference k-mers for each sample. The k-mer set of the reference needs to be computed only once.

`dicey chop -a -f hs37d5 -l 61 hs37d5.fa.gz`

```
dicey chop -a -f ${ID} -l 61 ${ID}.unitigs.fa.gz
sort -k1,1n -k2,2n -m <(zcat hs37d5.kmer.gz) <(zcat hs37d5.kmer.gz) <(zcat ${ID}.hashes.gz) | uniq -u | gzip -c > ${ID}.nonref.kmer.gz
```

### Non-reference k-mer counts

An example script to generate a count table of non-reference k-mers is included in the [aggegrate subdirectory](https://github.com/tobiasrausch/kmerdbg/tree/main/aggregate). Based on such a count table `nonref.count.gz`, one can create a classical site-frequency spectra - a log-log plot of the frequency of observed non-reference k-mer counts. Plotting scripts are in the [kmerfreq subdirectory](https://github.com/tobiasrausch/kmerdbg/tree/main/kmerfreq). 

![Non-reference k-mers](https://github.com/tobiasrausch/kmerdbg/blob/main/kmerfreq/vac.png?raw=true)

Using the [1000 Genomes sample table](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel), one can extract non-reference k-mers by population or singleton k-mers and population-specific k-mers ([stats subdirectory](https://github.com/tobiasrausch/kmerdbg/tree/main/stats)).

![Population non-ref k-mers](https://github.com/tobiasrausch/kmerdbg/blob/main/stats/nonref.jpg?raw=true)

### QV assembly quality estimation

With the [unitig files](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.0/dbg_1kgp/) you can also estimate the assembly QV quality of matched samples using [yak](https://github.com/lh3/yak).

```
yak count -k 31 -K1.5g -t 16 -o sample.dbg.yak sample.unitigs.fa.gz
yak qv -t 16 -p -K 1.5g -l 100k sample.dbg.yak sample.assembly.h1.fa > hap1.yak.qv.seq.txt
yak qv -t 16 -p -K 1.5g -l 100k sample.dbg.yak sample.assembly.h2.fa > hap2.yak.qv.seq.txt
```

Citation
--------

The compacted de Bruijn graphs have been built as part of the 1000 Genomes long-read sequencing project. If you use them please cite

Siegfried Schloissnig, Samarendra Pani, Bernardo Rodriguez-Martin, Jana Ebler, Carsten Hain, Vasiliki Tsapalou, Arda Söylev, Patrick Hüther, Hufsah Ashraf, Timofey Prodanov, Mila Asparuhova, Sarah Hunt, Tobias Rausch, Tobias Marschall, Jan O Korbel.              
Long-read sequencing and structural variant characterization in 1,019 samples from the 1000 Genomes Project.
bioRxiv. 2024 Apr 20:2024.04.18.590093.
[https://doi.org/10.1101/2024.04.18.590093](https://doi.org/10.1101/2024.04.18.590093)

Credits
-------
[HTSlib](https://github.com/samtools/htslib) and [samtools](https://github.com/samtools/samtools) for genomic data processing, [Boost](https://www.boost.org/) for various data structures and algorithms, [Lighter](https://github.com/mourisl/Lighter) for sequencing error correction and [BCALM 2](https://github.com/GATB/bcalm) for the compacted de Bruijn graph construction.
