''' This file downloads all nessesary data from the internet to ../../data. '''
'''
	2014-04-13
	* Update to accomodate different species and ensembl releases
	* Remove non chromosome entries from GTF
	* The script works for any organism with [0-9]*, X, Y, W, Z, MT chromosomes
	Execution:
		python customize_download.py -o homo_sapiens -r 75

	Last Update: Chelsea Ju
'''

import os, argparse, re, datetime
import generate_fasta

def echo(message):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(message))

def download(release, organism):

	# specify directory for givenc release and organism
	if(release == "current"):
		gtf_dir = "ftp://ftp.ensembl.org/pub/current_gtf/" + organism
		gvf_dir = "ftp://ftp.ensembl.org/pub/current_variation/gvf/" + organism
		fasta_dir = "ftp://ftp.ensembl.org/pub/current_fasta/" + organism + "/dna"
	else:
		gtf_dir = "ftp://ftp.ensembl.org/pub/release-%s/gtf/%s" %(release, organism)
		gvf_dir = "ftp://ftp.ensembl.org/pub/release-%s/variation/gvf/%s" %(release, organism)
		fasta_dir = "ftp://ftp.ensembl.org/pub/release-%s/fasta/%s/dna/" %(release, organism)


	root = "../../data/" + organism + "/" + release + "/"
        os.system("rm -rf " + root)
        os.system("mkdir -p " + root)
	os.chdir(root)

	## GTF
	os.system("wget %s/*.gtf.gz" %(gtf_dir))

	# retreive the prefix of the file (based on the latest download filename)
	gtf_file = max([f for f in os.listdir('.') if f.lower().endswith('.gtf.gz')],
                       key=os.path.getctime)
	prefix_search = re.search('(.*).gtf', gtf_file)
	if(not prefix_search):
		print("Could not find the gtf.gz file!")
	 	sys.exit(0)

	prefix = prefix_search.group(1)
	# unzip GTF file
	os.system("gunzip %s " % (prefix + ".gtf.gz"))

	# extract protein coding entries
        gtf_origin = prefix + ".gtf"
        gtf_protein_only_origin = prefix+".protein_coding.gtf"
        os.system("grep 'protein_coding' %s > %s" % (gtf_origin, gtf_protein_only_origin + ".backup"))

        gtf_backup = gtf_origin + ".backup"
	# remove none chromosome entries (not in 1-19, X, Y, MT) from the GTF
	os.system("grep -E '^MT|^[XY]|^[0-9]*\t' %s > %s" %
                  (gtf_protein_only_origin + ".backup", gtf_protein_only_origin))
	os.system("mv %s %s" % (gtf_origin, gtf_origin + ".backup"))
	os.system("grep -E '^MT|^[XY]|^[0-9]*\t' %s > %s" %(gtf_origin + ".backup", gtf_origin))

	## GVF
	gvf_file = organism.capitalize() + ".gvf.gz"
	# os.system("wget %s/%s" %(gvf_dir, gvf_file))
	# os.system("gunzip %s " %(gvf_file))

	## DNA sequence
	target = prefix + ".dna.fa"
	os.system("wget %s/%s" %(fasta_dir, prefix+".dna.chromosome.[0-9]*.fa.gz"))
	os.system("wget %s/%s" %(fasta_dir, prefix+".dna.chromosome.[XYWZ].fa.gz"))
	os.system("wget %s/%s" %(fasta_dir, prefix+".dna.chromosome.MT.fa.gz"))
	os.system("gunzip %s" %(prefix + ".dna.chromosome.*.fa.gz"))
	os.system("cat %s >> %s" %(prefix + ".dna.chromosome.*.fa", target))
	os.system("mkdir -p fa")
	os.system("mv %s fa/" %(prefix + ".dna.chromosome.*.fa"))
	os.chdir("fa")
	# change filenames
	fa_file = [f for f in os.listdir('.') if f.lower().endswith('.fa')]
	for f in fa_file:
		fa_prefix = re.search('.*.chromosome.(.*).fa', f)
		if(fa_prefix):
			print("mv %s %s" %(f, fa_prefix.group(1)+".fa"))
			os.system("mv %s %s" %(f, fa_prefix.group(1)+".fa"))
	os.chdir("..")
        return root, gtf_origin, gtf_protein_only_origin, gvf_file,  target


def main(parser):

    options = parser.parse_args()
    release = options.release
    organism = options.organism

    echo("Download Data For %s release %s" %(organism, release))
    root, gtf_file, gtf_protein_coding_file, gvf_file, fa_file = download(release, organism)
    echo("Donwload Completed")
    transcript_fasta, gene_fasta, transcript_fasta_pc, gene_fasta_pc = generate_fasta.generate(
        root, gtf_file, gtf_protein_coding_file, fa_file)
    print transcript_fasta, gene_fasta, transcript_fasta_pc, gene_fasta_pc

if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='customize_dowload.py')
    parser.add_argument("-r", "--release", dest="release", type=str,
                        help="release version [ex 75], default = current release",
                        required = False, default="current")
    parser.add_argument("-o", "--organism", dest="organism", type=str,
                        help="species name [ex homo_sapiens], default = mus_musculus (mouse)",
                        required = False, default = "mus_musculus")

    main(parser)
