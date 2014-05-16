RNA-Skim
========

RNA-Skim: a rapid method for RNA-Seq quantification at transcript level


How to compile RNA-Skim?
------------------------

RNA-Skim is implemented in C++ (heavily using C++11 standard). Please make sure that g++ (>= 4.7) is installed. Note: The default compiler of MacOS is clang, and currently RNA-Skim cannot be compiled by clang, so, please make sure that g++ is your default compiler, e.g., "export CXX=/opt/local/bin/g++-mp-4.8". If you set the default compiler to g++, please run the following commands to compile the executables:


```bash
./prepare.sh
cd src
make all
```

You can run "find *_test  -exec ./{} \;" in the src folder to test RNA-Skim, and if all tests are passed, you have successfully compiled RNA-Skim.

Workflow of RNA-Skim
--------------------

The workflow of RNA-Skim includes two parts: the preparation stage (time-consuming; excutables include: rs\_cluster, rs\_index, rs\_select), and the quantification stage (rapid; excutables include: rs\_count and rs\_estimate). We will go through the usage of each executable in the following. 

You can also use "--helpsort" parameter to check all parameters supported by these executables, e.g.,

```bash
./rs_cluster --helpsort
```

Download data from Ensembl
--------------------------

We have a python script to help users to download data from Ensembl. Please go to the src/prepare folder, and run the following commands:

```
cd src/prepare
python customize_download.py -o homo_sapiens -r current
python customize_download.py -o mus_musculus -r current
```

And the "-o" parameter means the population type, and the "-r" means the release number. The data will be downloaded at "data/homo\_sapiens/current" and "data/mus\_musculus/current". At the end, there is a simple description of how to prepare the "clustered\_gene.fa.sk" for quantification, and how to analyze a given set of fasta files. You can also take a look at the following document to check the usage of each command.

rs_cluster
----------

rs_cluster clusters the similar genes based on their sequence similarity. 

Before this workflow, you need to prepare a specialized FASTA format file, which contains the transcriptome of the subject of interest. If you use customize\_download.py, at the end, it outputs two specialized FASTA format files: one contains all transcripts, and the other only contains the protein coding transcripts (the latter is recommended). 

This part explains what is the difference between the regular FASTA format and the specialized FASTA format. Assuming we have two genes G1 and G2, and G1 has two transcripts T1 and T2, and G2 has one transcript T3. And the sequences of these transcripts are: T1: ATTA, T2: GAGA, T3: TTAA.

You need to prepare a fasta file looks like this:
```
>G1|T1|T2
ATTA|GAGA
>G2|T3
TTAA
```
The id lines start with the gene id, followed by the transcript ids, seperated by "|". The sequence lines contain the sequences of transcripts, and are also seperated by "|".

If you want to do the transcript-based clustering (Not recommended), the fasta file looks like this:
```
>T1|T1
ATTA
>T2|T2
GAGA
>T3|T3
TTAA
```

There shuold be no duplications on the first field. 

Let us assume the fasta file is named "gene.fa", and you can use the following command to generate the clustering result. 

```bash
GLOG_logtostderr=1 ./rs_cluster  -gene_fasta=gene.fa -num_threads=4 -output=clustered.fa -rs_length=60
```

The rs_length parameter is the length of k-mer used for calculating the similarity. 
And the clustered.fa is also in the same specialized FASTA format. In this case, each item represents a cluster, and the first field is randomly selected from the genes in the cluster (we do not track the gene id in the future analysis). 

rs_index
--------

Now you can run the following command to find all sig-mer regions.

```
GLOG_logtostderr=1  ./rs_index -transcript_fasta=clustered.fa -index_file=clustered_gene.fa.pb -rs_length=60 -num_threads 4
```
A sig-mer region is a sequence that all k-mers from the region are sig-mers of the transcript cluster. (Check out the document for the GeneSignatures class at rnasigs.proto for more details). The clustered\_gene.fa.pb contains the corresponding GeneSignatures class for every transcript cluster. Please make sure to use the same value for rs\_length for all executables. 

rs_select
---------

OK, now you have the clustered_gene.fa.pb file. Let's select sig-mers from all sig-mer regions:

```GLOG_logtostderr=1 ./rs_select -index_file=clustered_gene.fa.pb -selected_keys_file=clustered_gene.fa.sk```

**If you use sig-mer size other than the default value, you should specify the length in the parameter list as well,** for example:

```GLOG_logtostderr=1 ./rs_select -index_file=clustered_gene.fa.pb -selected_keys_file=clustered_gene.fa.sk  -rs_length=60```


You may see some warnings like this,

```E0309 21:20:25.820291 1990034192 rs_select.cc:130] ENSMUST00000114890 has 1 rna_signatures. Skipped.```

For some transcripts, RNA-Skim cannot get enough number of sig-mers, so it does not quantify such transcripts because the results of those transcripts are not reliable. Most of such transcripts are either too short (less than 150 bps) or categorized as "predicted gene".

The output file (clustered_gene.fa.sk) is a list of SelectedKey objects. For details, please checkout the comments for SelectedKey at rnasigs.proto.

rs_count
--------

rs_count counts the occurrences of the sig-mers for a given RNA-Seq dataset

Now, we run rs\_count to count all sig-mers stored in the clustered_gene.fa.sk file.

```
GLOG_logtostderr=1  ../src/rs_count  -selected_keys_file=clustered_gene.fa.sk -count_file=clustered_gene.fa.cf -read_files1=../test/test.fastq_1 -read_files2=../test/test.fastq_2 -num_threads=1
```

This generates clustered\_gene.fa.cf file, which is almost identical with the clustered\_gene.fa.sk file, but the count fields in the SelectedKey object in the clustered_gene.fa.cf file is the real occurrences of their corresponding sig-mers.

rs_estimate
-----------

rs_estimate quantifies the abundances of transcripts based on the occurrences of sig-mers.

This command quantifies the transcriptome based on the counts of sig-mers in the clustered_gene.fa.cf file.
```
../src/rs_estimate -count_file=clustered_gene.fa.cf > estimation
```

There are four columns in the estimation file: transcript id; the length of the transcript; estimated number of relative reads; RPKM value of the transcript.

