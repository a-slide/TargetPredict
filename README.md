# TargetPredict 0.1

**Find the targets of a short RNA in the transcriptome using BLASTN and MIRANDA**

---

**Creation : 2015/03/08**

**Last update : 2015/03/25**

---

## Motivation and Principle

This program was develop to detect imperfect off-target effects of a small RNA in the human transcriptome.

#### GffFastaExtractor.py

Extract the sequences from a fasta file corresponding to features location indicated in a gff file. Features in the gff files can be restricted to a list of chromosomes and a list of types. Features can be extended in 5' and 3' and overlapping features can be merged together.

#### TargetPredict.py

Search targets for a query fasta sequence (short RNA) in a subject fasta sequence using Blastn and Miranda algorithms. The ouput of the program is parsed and returned as csv tables sorted according to the alignment score. 

## Dependencies

The program was developed under Linux Mint 17.2 and was not tested with other OS. In addition to a python2.7 envronment the following dependencies are required for proper program execution:

**GffFastaExtractor.py**

* Python package:  Biopython 1.65 +

**TargetPredict.py**

* NCBI Blast 2.2.28+
* Miranda v3.3a+

## Get and install TargetPredict

* Clone the repository ```git clone https://github.com/a-slide/TargetPredict.git```

* Enter the src folder of the program folder

* Make the main scripts executable ```sudo chmod u+x TargetPredict.py GffFastaExtractor.py```

* Finally, add or link TargetPredict.py and GffFastaExtractor.py to your PATH

## Usage

**GffFastaExtractor.py**

```
Usage: GffFastaExtractor.py -f FASTA_file -g GFF_file [...]

Options:
  --version             	show program's version number and exit
  -h, --help            	show this help message and exit
  -f FASTA, --fasta=FASTA	Path to the fasta file contaning the complete genome (can be gzipped)
  -g GFF, --gff=GFF     	Path to the gff file containing annotations of the genome file (can be gzipped)
  -o OFFSET, --offset=OFFSET	Bases to extract before and after the feature (default 0)
  --fusion              	Fuse overlapping features in a meta-feature (default False)
  --output_gff          	Output the gff file corresponding to the extracted features sequences (default False)
  --features=FEATURES   	Restrict extraction to a list of features. The list must be SPACE separated and quoted (default exon)
  --chromosomes=CHROMOSOMES 	Restrict extraction to a list of chromosomes. The list must be SPACE separated and quoted (default all)
```

**TargetPredict.py**

```
Usage: TargetPredict.py -s SUBJECT -q QUERY [...]

Options:
  --version             	show program's version number and exit
  -h, --help            	show this help message and exit
  -s SUBJECT, --subject=SUBJECT	Path to the subject in fasta
  -q QUERY, --query=QUERY	Path to the query in fasta
  --blastn_opt=BLASTN_OPT	Options to run blastn (default='-task blastn-short -evalue 100000000 -word_size 5 -max_target_seqs 10000')
  --miranda_opt=MIRANDA_OPT	Options to run miranda (default='-sc 160')
```

## Authors and Contact

Adrien Leger - 2015

* <adrien.leger@gmail.com> - <adrien.leger@inserm.fr> - <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089](http://www.atlantic-gene-therapies.fr/)
* [Drowned in genomics](http://a-slide.github.io/)
