#process: laminB1Lads.txt.gz model-based-cpg-islands-hg19.txt

all : laminB1Lads.txt.gz model-based-cpg-islands-hg19.txt
.PHONY : all

laminB1Lads.txt.gz:
	bash download.sh 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz'

model-based-cpg-islands-hg19.txt:
	bash download.sh 'http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt'

#cpg_islands_sorted.txt:
#	touch $@

#.INTERMEDIATE: cpg_islands_sorted.txt
