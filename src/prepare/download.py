''' This file downloads all nessesary data from the internet to ../../data. '''
import os


root = "../../data/"

os.system("mkdir -p " + root)

os.chdir(root)
os.system("wget ftp://ftp.ensembl.org/pub/release-70/gtf/mus_musculus/Mus_musculus.GRCm38.70.gtf.gz")
os.system("gunzip Mus_musculus.GRCm38.70.gtf.gz ")
os.system("grep 'protein_coding' Mus_musculus.GRCm38.70.gtf >  Mus_musculus.GRCm38.70.protein_conding.gtf")

os.system("wget ftp://ftp.ensembl.org/pub/release-71/variation/gvf/mus_musculus/Mus_musculus.gvf.gz")
os.system("gunzip Mus_musculus.gvf.gz ")

target = "Mus_musculus.GRCm38.70.dna.fa"
os.system("rm -f " + target)
for c in range(1, 20) + ['X', 'Y', 'MT']:
    ftp = "ftp://ftp.ensembl.org/pub/release-70/fasta/mus_musculus/dna/"
    f = "Mus_musculus.GRCm38.70.dna.chromosome." + str (c) + ".fa.gz"
    os.system("wget " + ftp + f)
    os.system("gunzip " + f)

os.system("cat *.GRCm38.70.dna.chromosome* >> " + target)

os.system("mkdir fa")
os.system("mv *chromosome* fa")
os.chdir("fa")

for c in range(1, 20) + ['X', 'Y', 'MT']:
    f = "Mus_musculus.GRCm38.70.dna.chromosome." + str (c) + ".fa"
    nf = str (c) + ".fa"
    os.system("mv " + f + " " + nf)

## TODO: remove none chromosome entries (not in 1-19, X, Y, MT) from the
## gtf file.
