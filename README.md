RNA-Skim
========

RNA-Skim: a rapid method for RNA-Seq quantification at transcript level


How to compile RNA-Skim?
------------------------

RNA-Skim is implemented in C++ (heavily using C++11 standard). Please make sure that g++ (>= 4.7) is installed. Note: The default compiler of MacOS is clang, and currently RNA-Skim cannot be compiled by clang, so, please make sure that g++ is your default compiler, e.g., "export CXX=/opt/local/bin/g++-mp-4.8". 



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


rs_cluster
----------

rs_cluster clusters the similar genes based on their sequence similarity. 

Before this workflow, you need to prepare a modified FASTA format file, which contains the transcriptome of the subject of interest. 

Assuming we have two genes G1 and G2, and G1 has two transcripts T1 and T2, and G2 has one transcript T3. And the sequences of these transcripts are: T1: ATTA, T2: GAGA, T3: TTAA.

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

