## JigsawSeq ##

Comprehensive sequence profiling of gene variant libraries is important as it is broadly applicable for fine mapping studies in genotype-phenotype relations. JigsawSeq is a protocol to generate and evaluate a full-length gene variant library by combining codon-based molecular barcoding strategy.
This program **JigsawSeq** is a package of perl scripts to perform de novo assembly with next-generation sequencing data from JigsawSeq experiments, and analyze the contigs.


### Getting Started ###

JigsawSeq package is consisted of stand-alone perl scripts and can be downloaded from [github](http://github.com/jy2/JigsawSeq)

    git clone https://github.com/jy2/JigsawSeq.git

Then, you can run **JigsawSeq** on cloned directory without installation.

Requirements:

- `perl` 5.10 or higher 
- `bwa` 0.79a
- `samtools` 0.1.19

cf. `bwa` and `samtools` are included in this package.  


Simply, **JigsawSeq** can be excuted like:

    ./Jigsaw.pl -F Ex1_f.fastq -R Ex1_r.fastq -V pBR322_vector.fasta -L 816 -O Ex1


The example fastq files (ex1\_f.fastq, ex1\_r.fastq) can be downloaded in here; [Ex1_f.fastq.gz](http://chem.yonsei.ac.kr/~duhee/files/Ex1_f.fastq.gz) and  [Ex1_r.fastq.gz](http://chem.yonsei.ac.kr/~duhee/files/Ex1_r.fastq.gz).

With default parameters, **JigsawSeq** will use ~20Gb of memory for analyzing a 1.5Gb size of fastq file (~5,000,000 reads for *KanR* library). But, the memory usage depends on several factors; amount of reads, length of contigs, mutation rate, and sequencing quality. Therefore, the exact memory usage will be varied at specific situations. We highly recommend to use Hiseq data instead of MiSeq or NextSeq data, since sequencing error rates of MiSeq/NextSeq data are generally higher than one of HiSeq data, and much higher memory space is required when using other machines. Note that you can reduce the memory usage by lowering `--bin_size` (default: 15,000,000) when running `Jigsaw.pl`.


The final analyzed results will be stored at ex1.contigs.result.This result file is a tab-delimited text file with columns; *contig name*, *length of contig*, *depth of contig*, *initial seed*, *terminal seed*, and *contig sequence*. The numbers on 4th and 5th column are matched to the seed sqeuences at the header. 

### Workflow of JigsawSeq ###

![WorkFlow](https://lh3.googleusercontent.com/-8jyrnSoEBt4/VTZfn4vR4-I/AAAAAAAAA3M/-tnhmpCXeiU/s512/WorkFlow.png)

### Commend-line options ###
#### Jigsaw.pl ####
- Description: main module to preprocessing raw fastq files and construct contigs. 
- Usage: Jigsaw.pl [OPTIONS] -F [input:fastq] -R [input:fastq] -V [input:fasta] -L [input:integer] -O [output:prefix]
- Parameters: Among the parameters below, -F, -R, -V, -O, -L are required parameters

    `-F, --forward` Input filename of forward reads (fastq format) 

    `-R, --reverse`   Input filename of reverse reads (fastq format) 

    `-V, --vector`    Input filename of vector sequence (fasta format). The representation of vector sequence must be followed; 3'side of the insertion site must be the first position of vector sequence and 5'-end of the insertion site must be the last position. 

    `-O, --output`    Prefix of temporary and final output files

    `-L, --length`    Expected length of contig will determine the maximum number of steps on exploring graph.

    `-k, --kmer`      Length of k-mer. Default = 120

    `-s, --step`    Step size (1-, 2-, or 3-mer) for exploring the graph. Default = 3

    `-m, --min_depth`   Minimum depth of nodes and edges. The nodes and edges with lower depths than `min_depth` will be excluded. Default = 2.

    `--cut_edge`   Cutoff ratio for edges. To save the exploring time and memory usage, and to avoid false positive contigs by errorous edges  the edges will be discarded if the depth of edges divided by the maximum depth of edges linked to the same node is lower than 1/`cut_edge`. Default = 150.

    `--min_seed`   Minumum depth of seeds. The seeds with lower depth than `min_seed` will be ignored. Default = 100.

    `--cut_seed`   Cutoff ratio for seeds. To avoid false positive contigs by errorous intitial/terminal seeds, the seeds will be neglected if the depth of seeds divided by the maximum depth of seeds is lower than 1/`cut_seed`. Default = 100.

    `-e, --exclude`  File name of problematic sequences (fasta format). Any nodes having these sequences will be excluded. Default = exclusion.fa

    `--read_length` Read length of raw reads. Default = 150

    `-t, --thread` Number of threads. Only for detecting seeds by using `bwa`. Default = 3.

    `-b, --bin_size` Number of reads in a bin. Lowering `--bin_size` will reduce the memory usage, but increase the running time. Default = 15000000.

    `-C, --cut_CV` Cutoff for coefficient of variation. This will be used only if the contigs were re-aligned by raw data. (-a option) Default = 0.2163.

    `-a, --realign`  Allow to realign raw reads to contigs. This mode will reduce false positive contigs but decrease sensitivity.

    `-v, --verbose`  Verbose-mode On. Temporary files will not be deleted.





### Information ###
- Version: r3
- Last modified: Apr-22-2015
- Contact: Jung-Ki Yoon M.D. (dr.jkyoon at gmail.com)
 



