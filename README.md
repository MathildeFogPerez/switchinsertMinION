# Switch region PIPELINE (MinION) #

Copyright (C) 2020  Mathilde Foglierini Perez

email: mathilde.perez@irb.usi.ch

### SUMMARY ###

We have made available here a series of scripts to analyze the Switch region (IGH locus) to find potential DNA insertions. 
This pipeline is an updated version of https://bitbucket.org/mathildefog/switchminion/src/master/.

The scripts are primarily intended as reference for manuscript "Antibody diversification through non-VDJ insertions" rather than a stand-alone application.

The input of the pipeline are 2D reads coming from target amplicon of the switch region sequenced by Oxford Nanopore Technologies, quality-based filtered by Metrichor (https://metrichor.com/s/). Data can be found at SRA db: PRJNA638005 accession number.

These scripts were run on Linux machines.


### LICENSES ###

This code is distributed open source under the terms of the GNU Free Documention License.


### INSTALL ###

Before the pipeline can be run, the following software are required:

a) LAST v963 http://last.cbrc.jp/

b) Bedtools v2.26 http://bedtools.readthedocs.io/en/latest/index.html#

c) BEDOPS v2.4.20 https://bedops.readthedocs.io/en/latest/

d) Java JDK 12 https://www.oracle.com/technetwork/java/javase/downloads/index.html


### PIPELINE ###

First download all scripts and needed files in a folder.
Then we create a directory for each donor/barcode ($DONOR) and move in the 2D passed fastq file ($DONOR.fastq).
Make the last index file of the human genome hg 38.
All the following command lines can be run in a bash script.

        $ DONOR="KdlR"
        $ FOLDER="/$DONORPATH/" 
        $ SCRIPTSFOLDER="/PATHTOSCRIPTS/"
        $ genome="/PATHtohg38LastIndex/hg38"
        $ gencode="/PATHTOSCRIPTS/gencode.v29.geneName.annotation.sorted.bed"
        $ SPECIES="human"	

        $ #Remove the reads that are below 700 bp
        $ awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 700) {print header, seq, qheader, qseq}}' < $DONOR.fastq > $DONOR.mini700bp.fastq
        
		$#Transform to fasta format
        $ sed -n '1~4s/^@/>/p;2~4p' $DONOR.mini700bp.fastq> $DONOR.mini700bp.fasta

        $ file=$DONOR.mini700bp.fasta
        $ last-train -P4 $genome $file > $DONOR.par
        $ lastal -P4 -p $DONOR.par $genome $file | last-split -m1e-6  > $DONOR.maf
        $ #Make Switch coordinates map for insert and non insert reads
        $ echo $DONOR.maf
        $ java -jar $SCRIPTSFOLDER/MakeSwitchCoordinatesMap.jar $DONOR $file $DONOR.maf $SPECIES -50
        $ sort-bed $DONOR.switchCoord.bed>$DONOR.switchCoord.sorted.bed
        $ bedmap --echo --echo-map --delim '\t' $DONOR.switchCoord.sorted.bed $gencode > $DONOR.switchCoord.sorted.annotated.bed
        $ #finally write back the annotation info in the tsv file
        $ java -jar $SCRIPTSFOLDER/AddAnnotationInfo.jar $DONOR
        $ sort-bed "$DONOR"_insertList_afterLAST.bed > "$DONOR"_insertList_afterLAST_sorted.bed
        $ bedtools merge -i "$DONOR"_insertList_afterLAST_sorted.bed -c 4 -o sum > "$DONOR"_FINALinsertList_afterLAST.tsv
        $ bedmap --echo --echo-map-id-uniq --delim '\t' "$DONOR"_FINALinsertList_afterLAST.tsv $gencode > "$DONOR"_FINALinsertList_AfterLAST_Annotated.bed
        $ echo "Number of inserts below"
        $ wc -l "$DONOR"_FINALinsertList_AfterLAST_Annotated.bed
        $ echo "Number of reads containing inserts: "
        $ grep ">" "$DONOR"_selectedReadsAfterLAST.fasta | wc -l


    A bed file with the insert coordinates and a fasta file with the reads containing inserts will be produced at the end of the pipeline.