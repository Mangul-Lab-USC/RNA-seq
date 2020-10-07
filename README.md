# RNA-seq

<table>
<thead>
<tr class="header">
<th><strong>Tool</strong></th>
<th><strong>Year Published</strong></th>
<th><strong>Notable Features</strong></th>
<th><strong>Programming language</strong></th>
<th><strong>Package manager</strong></th>
<th><strong>Required expertise</strong></th>
<th><strong>Software</strong></th>
<th><p><strong>Type of URL</strong></p>
<p><strong>1. Web services designed to host source code</strong></p>
<p><strong>2. Others (e.g personal and/or university web services)</strong></p></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td colspan="8">
<p><strong>a. Data quality control</strong></p>
</td>

</tr>
<tr class="even">
<td>iSeqQC<a href="#ref1"><sup>1</sup></a></td>
<td>2020</td>
<td>Expression-based raw data QC tool that detects outliers</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/gkumar09/iSeqQC">https://github.com/gkumar09/iSeqQC</a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>qsmooth<a href="#ref2"><sup>2</sup></a></td>
<td>2018</td>
<td>Adaptive smooth quantile normalization</td>
<td>R</td>
<td>Bioconductor</td>
<td>++</td>
<td><a href="http://bioconductor.org/packages/release/bioc/html/qsmooth.html"><span class="underline">http://bioconductor.org/packages/release/bioc/html/qsmooth.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>FastQC<a href="#ref3"><sup>3</sup></a></td>
<td>2018</td>
<td>Raw data QC tool for for high throughput sequence data</td>
<td>Java</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://urldefense.com/v3/__https://github.com/s-andrews/fastqc/__;!!LIr3w8kk_Xxm!_WedrTyYaJKECbSrD_jkoEg5ONGOvtR9H9CpJl1oLyf-ca_73n3OM8LfadvwIPw$"><span class="underline">https://github.com/s-andrews/fastqc/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>QC3<a href="https://paperpile.com/c/VjSfVH/oXXz9"><sup>4</sup></a></td>
<td>2014</td>
<td>Raw data QC tool detecting batch effect and cross contamination</td>
<td>Perl, R</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/slzhao/QC3"><span class="underline">https://github.com/slzhao/QC3</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>kPAL<a href="https://paperpile.com/c/VjSfVH/Z32Ou"><sup>5</sup></a></td>
<td>2014</td>
<td>Alignment-free assessment raw data QC tool by analyzing k-mer frequencies</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/LUMC/kPAL"><span class="underline">https://github.com/LUMC/kPAL</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>HTQC<a href="https://paperpile.com/c/VjSfVH/ZdZOL"><sup>6</sup></a></td>
<td>2013</td>
<td>Raw data QC read assessment and filtration</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/htqc/"><span class="underline">https://sourceforge.net/projects/htqc/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Trimmomatic<a href="https://paperpile.com/c/VjSfVH/euzDZ"><sup>7</sup></a></td>
<td>2014</td>
<td>Trimming of reads and removal of adapters</td>
<td>Java</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://www.usadellab.org/cms/index.php?page=trimmomatic"><span class="underline">http://www.usadellab.org/cms/index.php?page=trimmomatic</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>Skewer<a href="https://paperpile.com/c/VjSfVH/AYNs8"><sup>8</sup></a></td>
<td>2014</td>
<td>Adapter trimming of reads</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://sourceforge.net/projects/skewer"><span class="underline">https://sourceforge.net/projects/skewer</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Flexbar<a href="https://paperpile.com/c/VjSfVH/ZxxWC"><sup>9</sup></a></td>
<td>2012</td>
<td>Trimming of reads and adaptor removal</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/seqan/flexbar"><span class="underline">https://github.com/seqan/flexbar</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>QuaCRS<a href="https://paperpile.com/c/VjSfVH/9pl3i"><sup>10</sup></a></td>
<td>2014</td>
<td>Post QC tool by performing meta-analyses on QC metrics across large numbers of samples.</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/kwkroll32/QuaCRS"><span class="underline">https://github.com/kwkroll32/QuaCRS</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>BlackOPs<a href="https://paperpile.com/c/VjSfVH/kl83M"><sup>11</sup></a></td>
<td>2013</td>
<td>Post QC tool that simulates experimental RNA-seq derived from the reference genome and aligns these sequences and outputs a blacklist of positions and alleles caused by mismapping</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/rnaseqvariantbl/"><span class="underline">https://sourceforge.net/projects/rnaseqvariantbl/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>RSeQC<a href="https://paperpile.com/c/VjSfVH/hF0OY"><sup>12</sup></a></td>
<td>2012</td>
<td><p>Post QC evaluation of different aspects of RNA-seq</p>
<p>experiments, such as sequence quality, GC bias, nucleotide composition bias, sequencing depth, strand specificity, coverage uniformity and read distribution over the genome structure.</p></td>
<td>Python, C</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://rseqc.sourceforge.net/"><span class="underline">http://rseqc.sourceforge.net/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>RNA-SeQC<a href="https://paperpile.com/c/VjSfVH/UXSG9"><sup>13</sup></a></td>
<td>2012</td>
<td>RNA-seq metrics for post- quality control and process optimization</td>
<td>Java</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://software.broadinstitute.org/cancer/cga/rna-seqc"><span class="underline">https://software.broadinstitute.org/cancer/cga/rna-seqc</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>Seqbias<a href="https://paperpile.com/c/VjSfVH/x9hen"><sup>14</sup></a></td>
<td>2012</td>
<td>Post QC tool using a graphical model to increase accuracy of de novo gene annotation, uniformity of read coverage, consistency of nucleotide frequencies and agreement with qRT-PCR</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://master.bioconductor.org/packages/devel/bioc/html/seqbias.html"><span class="underline">http://master.bioconductor.org/packages/devel/bioc/html/seqbias.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>SAMStat<a href="https://paperpile.com/c/VjSfVH/iAYlH"><sup>15</sup></a></td>
<td>2011</td>
<td>Post QC tool which plotsPost nucleotide overrepresentation and other statistics in mapped and unmapped reads in a html page</td>
<td>C</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://samstat.sourceforge.net"><span class="underline">http://samstat.sourceforge.net</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Samtools<a href="https://paperpile.com/c/VjSfVH/VYwJx"><sup>16</sup></a></td>
<td>2009</td>
<td>Post QC tool using generic alignment format for storing read alignments against reference sequences and to visualize the Binary/Alignment Map (BAM).</td>
<td>C, Perl</td>
<td>Anaconda</td>
<td>+</td>
<td><a href="https://github.com/samtools/samtools"><span class="underline">https://github.com/samtools/samtools</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td colspan="8"><strong>b. Read alignment</strong></td>

</tr>
<tr class="odd">
<td>deSALT<a href="https://paperpile.com/c/VjSfVH/W81bL"><sup>17</sup></a></td>
<td>2019</td>
<td>Long transcriptomic read alignment with de Bruijn graph-based index</td>
<td>C</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/ydLiu-HIT/deSALT"><span class="underline">https://github.com/ydLiu-HIT/deSALT</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Magic-BLAST<a href="https://paperpile.com/c/VjSfVH/rXBdj"><sup>18</sup></a></td>
<td>2018</td>
<td>Aligner for long and short reads through optimization of a spliced alignment score</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://ncbi.github.io/magicblast/"><span class="underline">https://ncbi.github.io/magicblast/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Minimap2<a href="https://paperpile.com/c/VjSfVH/c2pYX"><sup>19</sup></a></td>
<td>2018</td>
<td>Alignment using seed chain alignment procedure</td>
<td>C, Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/lh3/minimap2"><span class="underline">https://github.com/lh3/minimap2</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>DART<a href="https://paperpile.com/c/VjSfVH/ilhnw"><sup>20</sup></a></td>
<td>2018</td>
<td>Burrows-Wheeler Transform based aligner which adopts partitioning strategy to divide a read into two groups</td>
<td>C/C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/lh3/minimap2"><span class="underline">https://github.com/hsinnan75/DART</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>MMR<a href="https://paperpile.com/c/VjSfVH/3n9g"><sup>21</sup></a></td>
<td>2016</td>
<td>Resolves the mapping location of multi-mapping reads, optimising for locally smooth coverage.</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/ratschlab/mmr"><span class="underline">https://github.com/ratschlab/mmr</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>ContextMap 2<a href="https://paperpile.com/c/VjSfVH/rrOR6"><sup>22</sup></a></td>
<td>2015</td>
<td>Allows parallel mapping against several reference genomes</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.bio.ifi.lmu.de/ContextMap"><span class="underline">http://www.bio.ifi.lmu.de/ContextMap</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>HISAT<a href="https://paperpile.com/c/VjSfVH/X61GA"><sup>23</sup></a></td>
<td>2015</td>
<td>Aligning reads using an indexing scheme based on the Burrows-Wheeler transform and the Ferragina-Manzini (FM) index</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://www.bio.ifi.lmu.de/ContextMap"><span class="underline">http://www.ccb.jhu.edu/software/hisat/index.shtml</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>Segemehl<a href="https://paperpile.com/c/VjSfVH/cBerW"><sup>24</sup></a></td>
<td>2014</td>
<td>Multi-split mapping for circular RNA, trans-splicing, and fusion events in addition to performing splice alignment</td>
<td>C, C++, Perl, Python, Shell (Bash)</td>
<td>Anaconda</td>
<td>++</td>
<td>(<a href="http://www.bioinf.uni-leipzig.de/Software/segemehl/"><span class="underline">http://www.bioinf.uni-leipzig.de/Software/segemehl/</span></a>).</td>
<td>2</td>
</tr>
<tr class="odd">
<td>JAGuaR<a href="https://paperpile.com/c/VjSfVH/Nhgde"><sup>25</sup></a></td>
<td>2014</td>
<td>Uses a modified GTF (Gene Transfer Format) of known splice sites to build the complete sequence from all reads mapped to the transcript.</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://www.bcgsc.ca/resources/software/jaguar"><span class="underline">https://www.bcgsc.ca/resources/software/jaguar</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>CRAC<a href="https://paperpile.com/c/VjSfVH/oInun"><sup>26</sup></a></td>
<td>2013</td>
<td>Uses double K-mer indexing and profiling approach to map reads, predict SNPs, gene fusions, repeat borders.</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://crac.gforge.inria.fr/"><span class="underline">http://crac.gforge.inria.fr/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>STAR<a href="https://paperpile.com/c/VjSfVH/HcGrE"><sup>27</sup></a></td>
<td>2013</td>
<td>Aligns long reads against genome reference database</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/alexdobin/STAR"><span class="underline">https://github.com/alexdobin/STAR</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Subread<a href="https://paperpile.com/c/VjSfVH/75Mhc"><sup>28</sup></a></td>
<td>2013</td>
<td>Mapping reads to a reference genome using multi-seed strategy, called seed-and-vote</td>
<td>C, R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/Rsubread.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/Rsubread.html</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>TopHat2<a href="https://paperpile.com/c/VjSfVH/v4QLu"><sup>29</sup></a></td>
<td>2013</td>
<td>Alignment of transcriptomes in the presence of insertions, deletions and gene fusions</td>
<td>C++, Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://ccb.jhu.edu/software/tophat/index.shtml"><span class="underline">http://ccb.jhu.edu/software/tophat/index.shtml</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>OSA<a href="https://paperpile.com/c/VjSfVH/6Rnd9"><sup>30</sup></a></td>
<td>2012</td>
<td>K-mer profiling approach to map reads</td>
<td>C#</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.arrayserver.com/wiki/index.php?title=OSA"><span class="underline">http://www.arrayserver.com/wiki/index.php?title=OSA</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>PASSion<a href="https://paperpile.com/c/VjSfVH/sfG8k"><sup>31</sup></a></td>
<td>2012</td>
<td>Pattern growth pipeline for splice junction detection</td>
<td>C++, Perl, Shell (Bash)</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://trac.nbic.nl/passion/"><span class="underline">https://trac.nbic.nl/passion/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>RUM<a href="https://paperpile.com/c/VjSfVH/FosW8"><sup>32</sup></a></td>
<td>2011</td>
<td>Comparative analysis of RNA-seq alignment algorithms and the RNA-seq unified mapper</td>
<td>Perl, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.cbil.upenn.edu/RUM/"><span class="underline">http://www.cbil.upenn.edu/RUM/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>SOAPSplice<a href="https://paperpile.com/c/VjSfVH/cabdO"><sup>33</sup></a></td>
<td>2011</td>
<td>Ab initio detection of splice junctions</td>
<td>Perl</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://subread.sourceforge.net/"><span class="underline">http://soap.genomics.org.cn/soapsplice.html</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>MapSplice<a href="https://paperpile.com/c/VjSfVH/L2tcL"><sup>34</sup></a></td>
<td>2010</td>
<td>De novo detection of splice junctions</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://subread.sourceforge.net/"><span class="underline">https://github.com/LiuBioinfo/MapSplice</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>SpliceMap<a href="https://paperpile.com/c/VjSfVH/yCVDK"><sup>35</sup></a></td>
<td>2010</td>
<td>De novo detection of splice junctions and RNA-seq alignment</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://web.stanford.edu/group/wonglab/SpliceMap/"><span class="underline">http://web.stanford.edu/group/wonglab/SpliceMap/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>Supersplat<a href="https://paperpile.com/c/VjSfVH/kY2xi"><sup>36</sup></a></td>
<td>2010</td>
<td>De novo detection of splice junctions</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://mocklerlab.org/tools/1/manual"><span class="underline">http://mocklerlab.org/tools/1/manual</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>HMMSplicer<a href="https://paperpile.com/c/VjSfVH/JrZZK"><sup>37</sup></a></td>
<td>2010</td>
<td>Detection of splice junctions of short sequence reads</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://derisilab.ucsf.edu/software/hmmsplicer"><span class="underline">http://derisilab.ucsf.edu/software/hmmsplicer</span></a></td>
<td><span class="underline">2</span></td>
</tr>
<tr class="even">
<td>QPALMA<a href="https://paperpile.com/c/VjSfVH/hGOfW"><sup>38</sup></a></td>
<td>2008</td>
<td>Spliced alignments of short sequence reads.</td>
<td>C++, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.raetschlab.org/suppl/qpalma"><span class="underline">http://www.raetschlab.org/suppl/qpalma</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td colspan="8"><strong>c. Gene annotations</strong></td>

</tr>
<tr class="even">
<td>SQANTI<a href="https://paperpile.com/c/VjSfVH/T0s0"><sup>39</sup></a></td>
<td>2018</td>
<td>Analyses quality of long reads transcriptomes and removes artefacts.</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/ConesaLab/SQANTI"><span class="underline">https://github.com/ConesaLab/SQANTI</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Annocript<a href="https://paperpile.com/c/VjSfVH/sAdu3"><sup>40</sup></a></td>
<td>2015</td>
<td>Databases are downloaded to annotate protein coding transcripts with the prediction of putative long non-coding RNAs in whole transcriptomes.</td>
<td>Perl, Python, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/frankMusacchia/Annocript"><span class="underline">https://github.com/frankMusacchia/Annocript</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>CIRI<a href="https://paperpile.com/c/VjSfVH/jtYuD"><sup>41</sup></a></td>
<td>2015</td>
<td>De novo circular RNA identification</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/ciri/"><span class="underline">https://sourceforge.net/projects/ciri/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>TSSAR<a href="https://paperpile.com/c/VjSfVH/16V6S"><sup>42</sup></a></td>
<td>2014</td>
<td>Automated de novo TSS annotation from differential RNA-seq data</td>
<td>Java, Perl, R</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://rna.tbi.univie.ac.at/TSSAR"><span class="underline">http://rna.tbi.univie.ac.at/TSSAR</span></a></td>
<td><span class="underline">2</span></td>
</tr>
<tr class="even">
<td><strong>d. Transcriptome assembly</strong></td>

</tr>
<tr class="odd">
<td>FLAIR<a href="https://paperpile.com/c/VjSfVH/TTIUT"><sup>43</sup></a></td>
<td>2020</td>
<td>Full-length alternative isoform analysis of RNA</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/BrooksLabUCSC/FLAIR">https://github.com/BrooksLabUCSC/FLAIR</a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Scallop<a href="https://paperpile.com/c/VjSfVH/vP2QA"><sup>44</sup></a></td>
<td>2017</td>
<td>Splice-graph-decomposition algorithm which optimizes two competing objectives while satisfying all phasing constraints posed by reads spanning multiple vertices</td>
<td>C++</td>
<td>Anaconda</td>
<td>+++</td>
<td><a href="https://urldefense.com/v3/__https://github.com/Kingsford-Group/scallop__;!!LIr3w8kk_Xxm!-CgzZgIIDwrzc-GOnR89b3NzWh3bCOxrry6bhmysxWqBenbkuAT9iGeajMOMUSw$"><span class="underline">https://github.com/Kingsford-Group/scallop</span></a></td>
<td><span class="underline">1</span></td>
</tr>
<tr class="odd">
<td>CLASS2<a href="https://paperpile.com/c/VjSfVH/l0gw"><sup>45</sup></a></td>
<td>2016</td>
<td>Splice variant annotation</td>
<td>C++, Perl, Shell</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://sourceforge.net/projects/splicebox/"><span class="underline">https://sourceforge.net/projects/splicebox/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>StringTie<a href="https://paperpile.com/c/VjSfVH/hhaWZ"><sup>46</sup></a></td>
<td>2015</td>
<td>Applies a network flow algorithm originally developed in optimization theory, together with optional de novo assembly, to assemble transcripts</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://ccb.jhu.edu/software/stringtie"><span class="underline">http://ccb.jhu.edu/software/stringtie</span></a></td>
<td><span class="underline">2</span></td>
</tr>
<tr class="odd">
<td>Bridger<a href="https://paperpile.com/c/VjSfVH/dTtcj"><sup>47</sup></a></td>
<td>2015</td>
<td>De novo transcript assembler using a mathematical model, called the minimum path cover</td>
<td>C++, Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/rnaseqassembly/files/?source=navbar"><span class="underline">https://sourceforge.net/projects/rnaseqassembly/files/?source=navbar</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Bayesembler<a href="https://paperpile.com/c/VjSfVH/qvR6K"><sup>48</sup></a></td>
<td>2014</td>
<td>Reference genome guided transcriptome assembly built on a Bayesian model</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/bioinformatics-centre/bayesembler"><span class="underline">https://github.com/bioinformatics-centre/bayesembler</span></a>.</td>
<td><span class="underline">1</span></td>
</tr>
<tr class="odd">
<td>SEECER<a href="https://paperpile.com/c/VjSfVH/Ql98S"><sup>49</sup></a></td>
<td>2013</td>
<td>De novo transcriptome assembly using hidden Markov Model (HMM) based method</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://sb.cs.cmu.edu/seecer/"><span class="underline">http://sb.cs.cmu.edu/seecer/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>BRANCH<a href="https://paperpile.com/c/VjSfVH/BY7aZ"><sup>50</sup></a></td>
<td>2013</td>
<td>De novo transcriptome assemblies by using genomic information that can be partial or complete genome sequences from the same or a related organism.</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/baoe/BRANCH"><span class="underline">https://github.com/baoe/BRANCH</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>EBARDenovo<a href="https://paperpile.com/c/VjSfVH/ZGu0v"><sup>51</sup></a></td>
<td>2013</td>
<td>De novo transcriptome assembly uses an efficient chimera-detection function</td>
<td>C#</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/ebardenovo/"><span class="underline">https://sourceforge.net/projects/ebardenovo/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Oases<a href="https://paperpile.com/c/VjSfVH/ZrrJT"><sup>52</sup></a></td>
<td>2012</td>
<td>De novo transcriptome assembly using k-mer profiling and building a de Brujin graph</td>
<td>C</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/dzerbino/oases/tree/master"><span class="underline">https://github.com/dzerbino/oases/tree/master</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Cufflinks<a href="https://paperpile.com/c/VjSfVH/bhdyI"><sup>53</sup></a></td>
<td>2012</td>
<td>Ab initio transcript assembly, estimates their abundances, and tests for differential expression</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/cole-trapnell-lab/cufflinks"><span class="underline">https://github.com/cole-trapnell-lab/cufflinks</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>IsoInfer<a href="https://paperpile.com/c/VjSfVH/S3Nqy"><sup>54</sup></a></td>
<td>2011</td>
<td>Infer isoforms from short reads</td>
<td>C/C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.cs.ucr.edu/~jianxing/IsoInfer.html"><span class="underline">http://www.cs.ucr.edu/~jianxing/IsoInfer.html</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>IsoLasso<a href="https://paperpile.com/c/VjSfVH/0Y01n"><sup>55</sup></a></td>
<td>2011</td>
<td>Reference genome guided using LASSO regression approach</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://alumni.cs.ucr.edu/~liw/isolasso.html"><span class="underline">http://alumni.cs.ucr.edu/~liw/isolasso.html</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>Trinity<a href="https://paperpile.com/c/VjSfVH/Klb2k"><sup>56</sup></a></td>
<td>2011</td>
<td>De novo transcriptome assembly</td>
<td>C++, Java, Perl, R, Shell (Bash)</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/trinityrnaseq/trinityrnaseq/wiki"><span class="underline">https://github.com/trinityrnaseq/trinityrnaseq/wiki</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Trans-ABySS<a href="https://paperpile.com/c/VjSfVH/Y2S4e"><sup>57</sup></a></td>
<td>2010</td>
<td>De novo short-read transcriptome assembly and can also be used for fusion detection</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/bcgsc/transabyss"><span class="underline">https://github.com/bcgsc/transabyss</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Scripture<a href="https://paperpile.com/c/VjSfVH/NbxHf"><sup>58</sup></a></td>
<td>2010</td>
<td>Ab initio reconstruction of transcriptomes of pluripotent and lineage committed cells</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.broadinstitute.org/software/Scripture/"><span class="underline">www.broadinstitute.org/software/Scripture/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td colspan="8"><strong>e. Transcriptome quantification</strong></td>

</tr>
<tr class="even">
<td>TALON<a href="https://paperpile.com/c/VjSfVH/Y1mRt"><sup>59</sup></a></td>
<td>2019</td>
<td>Long-read transcriptome discovery and quantification</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/dewyman/TALON"><span class="underline">https://github.com/dewyman/TALON</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Salmon<a href="https://paperpile.com/c/VjSfVH/xwlWH"><sup>60</sup></a></td>
<td>2017</td>
<td>Composed of: lightweight-mapping model, an online phase that estimates initial expression levels and model parameters, and an offline phase that refines expression estimates models, and mesures sequence-specific, fragment GC, and positional biases</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/COMBINE-lab/Salmon"><span class="underline">https://github.com/COMBINE-lab/Salmon</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Kallisto<a href="https://paperpile.com/c/VjSfVH/nUid2"><sup>61</sup></a></td>
<td>2016</td>
<td>K-mer based pseudoalignment for allligment free transcript and gene expression quantification</td>
<td>C, C++, Perl</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/pachterlab/kallisto"><span class="underline">https://github.com/pachterlab/kallisto</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Wub<a href="https://paperpile.com/c/VjSfVH/ZSFmM"><sup>62</sup></a></td>
<td>2016</td>
<td>Sequence and error simulation tool to calculate read and genome assembly accuracy.</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/nanoporetech/wub"><span class="underline">https://github.com/nanoporetech/wub</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Rcount<a href="https://paperpile.com/c/VjSfVH/gYga0"><sup>63</sup></a></td>
<td>2015</td>
<td>GUI based tool used for quantification using counts per feature</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://github.com/MWSchmid/Rcount"><span class="underline">https://github.com/MWSchmid/Rcount</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Ht-seq<a href="https://paperpile.com/c/VjSfVH/dNzvY"><sup>64</sup></a></td>
<td>2015</td>
<td>Calculates gene counts by counting number of reads overlapping genes</td>
<td>Python</td>
<td>pip</td>
<td>++</td>
<td><a href="https://htseq.readthedocs.io/en/release_0.11.1/overview.html"><span class="underline">https://htseq.readthedocs.io/en/release_0.11.1/overview.html</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>EMSAR<a href="https://paperpile.com/c/VjSfVH/PdqIQ"><sup>65</sup></a></td>
<td>2015</td>
<td>Estimation by mappability-based segmentation and reclustering using a joint Poisson model</td>
<td>C</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/parklab/emsar"><span class="underline">https://github.com/parklab/emsar</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Maxcounts<a href="https://paperpile.com/c/VjSfVH/kujz1"><sup>66</sup></a></td>
<td>2014</td>
<td>Quantify the expression assigned to an exon as the maximum of its per-base counts</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://urldefense.com/v3/__http://sysbiobig.dei.unipd.it/?q=Software*MAXCOUNTS__;Iw!!LIr3w8kk_Xxm!4MLTfE0Z7OacNJeXRmNGQeUZyIJFyAEl1RVU-fe_NYWYJAuOfXufQPMlJas-SFM$"><span class="underline">http://sysbiobig.dei.unipd.it/?q=Software#MAXCOUNTS</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>FIXSEQ<a href="https://paperpile.com/c/VjSfVH/w2lKf"><sup>67</sup></a></td>
<td>2014</td>
<td>A nonparametric and universal method for processing per-base sequencing read count data.</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://bitbucket.org/thashim/fixseq/src/master/"><span class="underline">https://bitbucket.org/thashim/fixseq/src/master/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Sailfish<a href="https://paperpile.com/c/VjSfVH/CfcCV"><sup>68</sup></a></td>
<td>2014</td>
<td>EM based quantification using statistical coupling between k-mers.</td>
<td>C, C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/kingsfordgroup/sailfish"><span class="underline">https://github.com/kingsfordgroup/sailfish</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Casper<a href="https://paperpile.com/c/VjSfVH/8K1xC"><sup>69</sup></a></td>
<td>2014</td>
<td>Bayesian modeling framework to quantify alternative splicing.</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://www.bioconductor.org/packages/release/bioc/html/casper.html"><span class="underline">http://www.bioconductor.org/packages/release/bioc/html/casper.html</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>MaLTA<a href="https://paperpile.com/c/VjSfVH/rqlp"><sup>70</sup></a></td>
<td>2014</td>
<td>Simultaneous transcriptome assembly and quantification from Ion Torrent RNA-Seq data.</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://alan.cs.gsu.edu/NGS/?q=malta"><span class="underline">http://alan.cs.gsu.edu/NGS/?q=malta</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>Featurecounts<a href="https://paperpile.com/c/VjSfVH/7Qa39"><sup>71</sup></a></td>
<td>2014</td>
<td>Read summarization program for counting reads generated.</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://subread.sourceforge.net/"><span class="underline">http://subread.sourceforge.net/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>MITIE<a href="https://paperpile.com/c/VjSfVH/IBZ5"><sup>72</sup></a></td>
<td>2013</td>
<td>Transcript reconstruction and assembly from RNA-Seq data using mixed integer optimisation.</td>
<td>MATLAB, C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/ratschlab/MiTie"><span class="underline">https://github.com/ratschlab/MiTie</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>iReckon<a href="https://paperpile.com/c/VjSfVH/e16fg"><sup>73</sup></a></td>
<td>2013</td>
<td>EM-based method to accurately estimate the abundances of known and novel isoforms.</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://compbio.cs.toronto.edu/ireckon/"><span class="underline">http://compbio.cs.toronto.edu/ireckon/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>eXpress<a href="https://paperpile.com/c/VjSfVH/rKytO"><sup>74</sup></a></td>
<td>2013</td>
<td>Online EM based algorithm for quantification which considers one read at a time.</td>
<td>C++, Shell (Bash)</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://pachterlab.github.io/eXpress/manual.html"><span class="underline">https://pachterlab.github.io/eXpress/manual.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>BitSeq<a href="https://paperpile.com/c/VjSfVH/q2362"><sup>75</sup></a></td>
<td>2012</td>
<td>Bayesian transcript expression quantification and differential expression.</td>
<td>C++, R</td>
<td>Anaconda,Bioconductor</td>
<td>++</td>
<td><a href="http://bitseq.github.io/"><span class="underline">http://bitseq.github.io/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>IQSeq<a href="https://paperpile.com/c/VjSfVH/v6vlI"><sup>76</sup></a></td>
<td>2012</td>
<td>Integrated isoform quantification analysis.</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://archive.gersteinlab.org/proj/rnaseq/IQSeq/"><span class="underline">http://archive.gersteinlab.org/proj/rnaseq/IQSeq/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>CEM<a href="https://paperpile.com/c/VjSfVH/xC6tA"><sup>77</sup></a></td>
<td>2012</td>
<td>Statistical framework for both transcriptome assembly and isoform expression level estimation.</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://alumni.cs.ucr.edu/~liw/cem.html"><span class="underline">http://alumni.cs.ucr.edu/~liw/cem.html</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>SAMMate<a href="https://paperpile.com/c/VjSfVH/WDgEC"><sup>78</sup></a></td>
<td>2011</td>
<td>Analysis of differential gene and isoform expression.</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://sammate.sourceforge.net/"><span class="underline">http://sammate.sourceforge.net/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Isoformex<a href="https://paperpile.com/c/VjSfVH/lQkln"><sup>79</sup></a></td>
<td>2011</td>
<td>Estimation method to estimate the expression levels of transcript isoforms.</td>
<td>N/A</td>
<td>N/A</td>
<td>N/A</td>
<td><a href="http://bioinformatics.wistar.upenn.edu/isoformex">http://bioinformatics.wistar.upenn.edu/isoformex</a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>IsoEM<a href="https://paperpile.com/c/VjSfVH/XJjxm"><sup>80</sup></a></td>
<td>2011</td>
<td>EM based method for inference of isoform and gene-specific expression levels</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://dna.engr.uconn.edu/software/IsoEM/"><span class="underline">http://dna.engr.uconn.edu/software/IsoEM/</span></a>.</td>
<td>2</td>
</tr>
<tr class="even">
<td>RSEM<a href="https://paperpile.com/c/VjSfVH/YTV4V"><sup>81</sup></a></td>
<td>2011</td>
<td>Ab initio EM based method for inference of isoform and gene-specific expression levels</td>
<td>C++, Perl, Python, R</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/deweylab/RSEM"><span class="underline">https://github.com/deweylab/RSEM</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>EDASeq<a href="https://paperpile.com/c/VjSfVH/nZEoj"><sup>82</sup></a></td>
<td>2011</td>
<td><h3 id="within-lane-gc-content-normalization-between-sample-normalization-visualization.">Within-lane GC-content normalization, between-sample normalization, visualization.</h3></td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/devel/bioc/html/EDASeq.html"><span class="underline">https://bioconductor.org/packages/devel/bioc/html/EDASeq.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>MMSEQ<a href="https://paperpile.com/c/VjSfVH/iRLTs"><sup>83</sup></a></td>
<td>2011</td>
<td>Haplotype and isoform specific expression estimation</td>
<td>C++, R, Ruby, Shell (Bash)</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/eturro/mmseq"><span class="underline">https://github.com/eturro/mmseq</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>MISO<a href="https://paperpile.com/c/VjSfVH/hDihL"><sup>84</sup></a></td>
<td>2010</td>
<td>Statistical model that estimates expression of alternatively spliced exons and isoforms</td>
<td>C, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://miso.readthedocs.io/en/fastmiso/#latest-version-from-github"><span class="underline">https://miso.readthedocs.io/en/fastmiso/#latest-version-from-github</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>SOLAS<a href="https://paperpile.com/c/VjSfVH/GJ3dA"><sup>85</sup></a></td>
<td>2010</td>
<td>Prediction of alternative isoforms from exon expression levels</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://cmb.molgen.mpg.de/2ndGenerationSequencing/Solas/"><span class="underline">http://cmb.molgen.mpg.de/2ndGenerationSequencing/Solas/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>Rseq<a href="https://paperpile.com/c/VjSfVH/nCaHC"><sup>86</sup></a></td>
<td>2009</td>
<td>Statistical inferences for isoform expression</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www-personal.umich.edu/~jianghui/rseq/#download"><span class="underline">http://www-personal.umich.edu/~jianghui/rseq/#download</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>rQuant<a href="https://paperpile.com/c/VjSfVH/wuq2V"><sup>87</sup></a></td>
<td>2009</td>
<td>Estimating density biases and considering the read coverages at each nucleotide independently using quadratic programming</td>
<td>Matlab, Shell (Bash), Javascript</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://galaxy.inf.ethz.ch/?tool_id=rquantweb&amp;version=2.2&amp;__identifer=3iuqb8nb3wf"><span class="underline">https://galaxy.inf.ethz.ch/?tool_id=rquantweb&amp;version=2.2&amp;__identifer=3iuqb8nb3wf</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>ERANGE<a href="https://paperpile.com/c/VjSfVH/O6Lej"><sup>88</sup></a></td>
<td>2008</td>
<td>Mapping and quantifying mammalian transcripts</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://woldlab.caltech.edu/rnaseq"><span class="underline">http://woldlab.caltech.edu/rnaseq</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td colspan="8"><strong>f. Differential expression</strong></td>

</tr>
<tr class="odd">
<td>Swish<a href="https://paperpile.com/c/VjSfVH/Tn9qR"><sup>89</sup></a></td>
<td>2019</td>
<td>Non-parametric model for differential expression analysis using inferential replicate counts</td>
<td>R</td>
<td>Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/fishpond.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/fishpond.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Yanagi<a href="https://paperpile.com/c/VjSfVH/CAWU"><sup>90</sup></a></td>
<td>2019</td>
<td>Transcriptome segment analysis</td>
<td>Python/C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/HCBravoLab/yanagi"><span class="underline">https://github.com/HCBravoLab/yanagi</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Whippet<a href="https://paperpile.com/c/VjSfVH/s2eM"><sup>91</sup></a></td>
<td>2018</td>
<td>Quantification of transcriptome structure and gene expression analysis using EM.</td>
<td>Julia</td>
<td>N/A</td>
<td>+++</td>
<td><h3 id="httpsgithub.comtimbitzwhippet.jl"><a href="https://urldefense.com/v3/__https://github.com/timbitz/Whippet.jl__;!!LIr3w8kk_Xxm!79cxwkS1fWXUUrx2H80qWxFR7igz5ztodLUb5CVzMfe1UBF69YEXgU7ErLh3YeI$"><span class="underline">https://github.com/timbitz/Whippet.jl</span></a></h3></td>
<td>1</td>
</tr>
<tr class="even">
<td>ReQTL<a href="https://paperpile.com/c/VjSfVH/RWWwc"><sup>92</sup></a></td>
<td>2018</td>
<td>Identifies correlations between SNVs and gene expression from RNA-seq data</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/HorvathLab/ReQTL"><span class="underline">https://github.com/HorvathLab/ReQTL</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>vast-tools<a href="https://paperpile.com/c/VjSfVH/mkVp"><sup>93</sup></a></td>
<td>2017</td>
<td>Profiling and comparing alternative splicing events in RNA-Seq data and for downstream analyses of alternative splicing.</td>
<td>R, Perl</td>
<td>N/A</td>
<td>+++</td>
<td><h3 id="httpsgithub.comvastgroupvast-tools"><a href="https://urldefense.com/v3/__https://github.com/vastgroup/vast-tools__;!!LIr3w8kk_Xxm!79cxwkS1fWXUUrx2H80qWxFR7igz5ztodLUb5CVzMfe1UBF69YEXgU7Efajq4L4$"><span class="underline">https://github.com/vastgroup/vast-tools</span></a></h3></td>
<td>1</td>
</tr>
<tr class="even">
<td>Ballgown<a href="https://paperpile.com/c/VjSfVH/fWFKP"><sup>94</sup></a></td>
<td>2015</td>
<td>Linear model–based differential expression analyses</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://github.com/alyssafrazee/ballgown"><span class="underline">https://github.com/alyssafrazee/ballgown</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Limma/Voom<a href="https://paperpile.com/c/VjSfVH/BuZHd"><sup>95</sup></a></td>
<td>2014</td>
<td>Linear model-based differential expression and differential splicing analyses</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/limma.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/limma.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>rMATS<a href="https://paperpile.com/c/VjSfVH/0run"><sup>96</sup></a></td>
<td>2014</td>
<td>Detect major differential alternative splicing types in RNA-seq data with replicates.</td>
<td>Python, C++</td>
<td>Anaconda</td>
<td>++</td>
<td><h3 id="httprnaseq-mats.sourceforge.netrmats3.2.5"><a href="https://urldefense.com/v3/__http://rnaseq-mats.sourceforge.net/rmats3.2.5/__;!!LIr3w8kk_Xxm!79cxwkS1fWXUUrx2H80qWxFR7igz5ztodLUb5CVzMfe1UBF69YEXgU7EG33qPzA$"><span class="underline">http://rnaseq-mats.sourceforge.net/rmats3.2.5/</span></a></h3></td>
<td>1</td>
</tr>
<tr class="odd">
<td>DESeq2<a href="https://paperpile.com/c/VjSfVH/PKb3X"><sup>97</sup></a></td>
<td>2014</td>
<td>Differential analysis of count data, using shrinkage estimation for dispersions and fold changes</td>
<td>R</td>
<td>Bioconductor, CRAN</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/DESeq2.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Corset<a href="https://paperpile.com/c/VjSfVH/jAg8X"><sup>98</sup></a></td>
<td>2014</td>
<td>Differential gene expression analysis for de novo assembled transcriptomes</td>
<td>C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/Oshlack/Corset/wiki"><span class="underline">https://github.com/Oshlack/Corset/wiki</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>BADGE<a href="https://paperpile.com/c/VjSfVH/FU1Ma"><sup>99</sup></a></td>
<td>2014</td>
<td>Bayesian model for accurate abundance quantification and differential analysis</td>
<td>Matlab</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.cbil.ece.vt.edu/software.htm"><span class="underline">http://www.cbil.ece.vt.edu/software.htm</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>compcodeR<a href="https://paperpile.com/c/VjSfVH/AsPdk"><sup>100</sup></a></td>
<td>2014</td>
<td>Benchmarking of differential expression analysis methods</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><span class="underline"><a href="https://urldefense.com/v3/__https://www.bioconductor.org/packages/compcodeR/__;!!LIr3w8kk_Xxm!6TILZhniJETsClNPMtzy-2gCUBnzcxnShT0euBD2UWYBwbx7MrkpXun-DaiXzTc$"> https://www.bioconductor.org/packages/compcodeR/</a></span></td>
<td>1</td>
</tr>
<tr class="odd">
<td>metaRNASeq<a href="https://paperpile.com/c/VjSfVH/qlyC3"><sup>101</sup></a></td>
<td>2014</td>
<td>Differential meta-analyses of RNA-seq data</td>
<td>R</td>
<td>Anaconda, CRAN</td>
<td>++</td>
<td><a href="http://cran.r-project.org/web/packages/metaRNASeq"><span class="underline">http://cran.r-project.org/web/packages/metaRNASeq</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Characteristic Direction<a href="https://paperpile.com/c/VjSfVH/DUkqa"><sup>102</sup></a></td>
<td>2014</td>
<td>Geometrical multivariate approach to identify differentially expressed genes</td>
<td>R, Python, MATLAB</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://www.maayanlab.net/CD"><span class="underline">http://www.maayanlab.net/CD</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>HTSFilter<a href="https://paperpile.com/c/VjSfVH/KkVO6"><sup>103</sup></a></td>
<td>2013</td>
<td>Filter-replicated high-throughput transcriptome sequencing data</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://www.bioconductor.org/packages/release/bioc/html/HTSFilter.html"><span class="underline">http://www.bioconductor.org/packages/release/bioc/html/HTSFilter.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>NPEBSeq<a href="https://paperpile.com/c/VjSfVH/tePCt"><sup>104</sup></a></td>
<td>2013</td>
<td>Nonparametric empirical bayesian-based procedure for differential expression analysis</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://bioinformatics.wistar.upenn.edu/NPEBseq"><span class="underline">http://bioinformatics.wistar.upenn.edu/NPEBseq</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>EBSeq<a href="https://paperpile.com/c/VjSfVH/DFbQ5"><sup>105</sup></a></td>
<td>2013</td>
<td>Identifying differentially expressed isoforms.</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://bioconductor.org/packages/release/bioc/html/EBSeq.html"><span class="underline">http://bioconductor.org/packages/release/bioc/html/EBSeq.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>sSeq<a href="https://paperpile.com/c/VjSfVH/2YaM9"><sup>106</sup></a></td>
<td>2013</td>
<td>Shrinkage estimation of dispersion in Negative Binomial models</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://bioconductor.org/packages/release/bioc/html/sSeq.html"><span class="underline">http://bioconductor.org/packages/release/bioc/html/sSeq.html</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Cuffdiff2<a href="https://paperpile.com/c/VjSfVH/orKVr"><sup>107</sup></a></td>
<td>2013</td>
<td>Differential analysis at transcript resolution</td>
<td>C++, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/"><span class="underline">http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>SAMseq<a href="https://paperpile.com/c/VjSfVH/q2TVc"><sup>108</sup></a></td>
<td>2013</td>
<td>Nonparametric method with resampling to account for the different sequencing depths</td>
<td>R</td>
<td>CRAN</td>
<td>++</td>
<td><a href="https://rdrr.io/cran/samr/man/SAMseq.html"><span class="underline">https://rdrr.io/cran/samr/man/SAMseq.html</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>DSGseq <a href="https://paperpile.com/c/VjSfVH/3Ks2"><sup>109</sup></a></td>
<td>2013</td>
<td>NB-statistic method that can detect differentially spliced genes between two groups of samples without using a prior knowledge on the annotation of alternative splicing.</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://bioinfo.au.tsinghua.edu.cn/software/DSGseq/"><span class="underline">http://bioinfo.au.tsinghua.edu.cn/software/DSGseq/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>NOISeq<a href="https://paperpile.com/c/VjSfVH/lZsbu"><sup>110</sup></a></td>
<td>2011</td>
<td>Uses a non-parametric approach for differential expression analysis and can work in absence of replicates</td>
<td>R</td>
<td>Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/NOISeq.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/NOISeq.html</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>EdgeR<a href="https://paperpile.com/c/VjSfVH/pnOy0"><sup>111</sup></a></td>
<td>2010</td>
<td>Examining differential expression of replicated count data and differential exon usage</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/edgeR.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>DEGseq<a href="https://paperpile.com/c/VjSfVH/Z0J7"><sup>112</sup></a></td>
<td>2010</td>
<td>Identify differentially expressed genes or isoforms for RNA-seq data from different samples.</td>
<td>R</td>
<td>Bioconductor</td>
<td>++</td>
<td><a href="http://bioconductor.org/packages/release/bioc/html/DEGseq.html"><span class="underline">http://bioconductor.org/packages/release/bioc/html/DEGseq.html</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td colspan="8"><strong>g. RNA splicing</strong></p></td>

</tr>
<tr class="even">
<td>LeafCutter<a href="https://paperpile.com/c/VjSfVH/munu8"><sup>113</sup></a></td>
<td>2018</td>
<td>Detects differential splicing and maps quantitative trait loci (sQTLs).</td>
<td>R</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/davidaknowles/leafcutter.git"><span class="underline">https://github.com/davidaknowles/leafcutter.git</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>MAJIQ-SPEL<a href="https://paperpile.com/c/VjSfVH/2Y6cD"><sup>114</sup></a></td>
<td>2018</td>
<td>Visualization, interpretation, and experimental validation of both classical and complex splicing variation and automated RT-PCR primer design.</td>
<td>C++, Python</td>
<td>N/A</td>
<td>+++</td>
<td><strong><a href="https://urldefense.com/v3/__https://galaxy.biociphers.org/galaxy/root?tool_id=majiq_spel__;!!LIr3w8kk_Xxm!9HDV79ohdTaNgi9fLGjyLb4eIakLRhb1CPCRFebtWc9ahY-6GQpBMKvS0LkeG0E$"><span class="underline">https://galaxy.biociphers.org/galaxy/root?tool_id=majiq_spel</span></a></strong></td>
<td>2</td>
</tr>
<tr class="even">
<td>MAJIQ<a href="https://paperpile.com/c/VjSfVH/b9zV"><sup>115</sup></a></td>
<td>2016</td>
<td>Web-tool that takes as input local splicing variations (LSVs) quantified from RNA-seq data and provides users with a visualization package (VOILA) and quantification of gene isoforms.</td>
<td>C++, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://majiq.biociphers.org/commercial.php"><span class="underline">https://majiq.biociphers.org/commercial.php</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>SplAdder<a href="https://paperpile.com/c/VjSfVH/zvhJm"><sup>116</sup></a></td>
<td>2016</td>
<td>Identification, quantification, and testing of alternative splicing events</td>
<td>Python</td>
<td>PyPI</td>
<td>+++</td>
<td><a href="http://github.com/ratschlab/spladder">http://github.com/ratschlab/spladder</a></td>
<td>1</td>
</tr>
<tr class="even">
<td>SplicePie<a href="https://paperpile.com/c/VjSfVH/9WudK"><sup>117</sup></a></td>
<td>2015</td>
<td>Detection of alternative, non-sequential and recursive splicing</td>
<td>Perl, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/pulyakhina/splicing_analysis_pipeline"><span class="underline">https://github.com/pulyakhina/splicing_analysis_pipeline</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>SUPPA<a href="https://paperpile.com/c/VjSfVH/5JQWU"><sup>118</sup></a></td>
<td>2015</td>
<td>Alternative splicing analysis</td>
<td>Python, R</td>
<td>Anaconda</td>
<td>++</td>
<td><h3 id="httpsgithub.comcomprnasuppa"><a href="https://urldefense.com/v3/__https://github.com/comprna/SUPPA__;!!LIr3w8kk_Xxm!79cxwkS1fWXUUrx2H80qWxFR7igz5ztodLUb5CVzMfe1UBF69YEXgU7EbRA4E0w$"><span class="underline">https://github.com/comprna/SUPPA</span></a></h3></td>
<td>1</td>
</tr>
<tr class="even">
<td>SNPlice<a href="https://paperpile.com/c/VjSfVH/5NuGs"><sup>119</sup></a></td>
<td>2015</td>
<td>Identifying variants that modulate Intron retention</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://code.google.com/p/snplice/"><span class="underline">https://code.google.com/p/snplice/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>IUTA<a href="https://paperpile.com/c/VjSfVH/bt7r8"><sup>120</sup></a></td>
<td>2014</td>
<td>Detecting differential isoform usage</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://www.niehs.nih.gov/research/resources/software/biostatistics/iuta/index.cfm"><span class="underline">http://www.niehs.nih.gov/research/resources/software/biostatistics/iuta/index.cfm</span></a>.</td>
<td>1</td>
</tr>
<tr class="even">
<td>SigFuge<a href="https://paperpile.com/c/VjSfVH/dc2Lc"><sup>121</sup></a></td>
<td>2014</td>
<td>Identifying genomic loci exhibiting differential transcription patterns</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://bioconductor.org/packages/release/bioc/html/SigFuge.html"><span class="underline">http://bioconductor.org/packages/release/bioc/html/SigFuge.html</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>FineSplice<a href="https://paperpile.com/c/VjSfVH/itVqO"><sup>122</sup></a></td>
<td>2014</td>
<td>Splice junction detection and quantification</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/p/finesplice/"><span class="underline">https://sourceforge.net/p/finesplice/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>PennSeq<a href="https://paperpile.com/c/VjSfVH/c1nBL"><sup>123</sup></a></td>
<td>2014</td>
<td>Statistical method that allows each isoform to have its own non-uniform read distribution</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://sourceforge.net/projects/pennseq"><span class="underline">http://sourceforge.net/projects/pennseq</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>FlipFlop<a href="https://paperpile.com/c/VjSfVH/xFcS0"><sup>124</sup></a></td>
<td>2014</td>
<td>RNA isoform identification and quantification with network flows</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/flipflop.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/flipflop.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>GESS<a href="https://paperpile.com/c/VjSfVH/bg83O"><sup>125</sup></a></td>
<td>2014</td>
<td>Graph-based exon-skipping scanner for de novo detection of skipping event sites</td>
<td>N/A</td>
<td>N/A</td>
<td>N/A</td>
<td><a href="http://jinlab.net/GESS_Web/"><span class="underline">http://jinlab.net/GESS_Web/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>spliceR<a href="https://paperpile.com/c/VjSfVH/ze0g0"><sup>126</sup></a></td>
<td>2013</td>
<td>Classification of alternative splicing and prediction of coding potential</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://www.bioconductor.org/packages/2.13/bioc/html/spliceR.html"><span class="underline">http://www.bioconductor.org/packages/2.13/bioc/html/spliceR.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>RNASeq-MATS<a href="https://paperpile.com/c/VjSfVH/gNSyJ"><sup>127</sup></a></td>
<td>2013</td>
<td>Detects and analyzes differential alternative splicing events</td>
<td>C, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://rnaseq-mats.sourceforge.net/">http://rnaseq-mats.sourceforge.net/</a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>SplicingCompass<a href="https://paperpile.com/c/VjSfVH/iR3p2"><sup>128</sup></a></td>
<td>2013</td>
<td>Differential splicing detection</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://www.ichip.de/software/SplicingCompass.html"><span class="underline">http://www.ichip.de/softwa</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>DiffSplice<a href="https://paperpile.com/c/VjSfVH/36q9H"><sup>129</sup></a></td>
<td>2013</td>
<td>Genome-wide detection of differential splicing</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.netlab.uky.edu/p/bioinfo/DiffSplice"><span class="underline">http://www.netlab.uky.edu/p/bioinfo/DiffSplice</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>DEXSeq<a href="https://paperpile.com/c/VjSfVH/TbBMF"><sup>130</sup></a></td>
<td>2012</td>
<td>Statistical method to test for differential exon usage.</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/DEXSeq.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/DEXSeq.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>SpliceSeq<a href="https://paperpile.com/c/VjSfVH/WUPzf"><sup>131</sup></a></td>
<td>2012</td>
<td>Identifies differential splicing events between test and control groups.</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://bioinformatics.mdanderson.org/main/SpliceSeq:Overview"><span class="underline">http://bioinformatics.mdanderson.org/main/SpliceSeq:Overview</span></a>.</td>
<td>2</td>
</tr>
<tr class="odd">
<td>JuncBASE<a href="https://paperpile.com/c/VjSfVH/q9K4"><sup>132</sup></a></td>
<td>2011</td>
<td>Identification and quantification of alternative splicing, including unannotated splicing</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://urldefense.com/v3/__https://github.com/anbrooks/juncBASE__;!!LIr3w8kk_Xxm!54q6C4EEhOKR6mGRcUzOweyB7vGk_4WfLJskEmh_sOV2eEfdSqMZ9G3WH65TQ2M$"><span class="underline">https://github.com/anbrooks/juncBASE</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>ALEXA-seq<a href="https://paperpile.com/c/VjSfVH/CvwSl"><sup>133</sup></a></td>
<td>2010</td>
<td>Alternative expression analysis.</td>
<td>Perl, R, Shell (Bash)</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.alexaplatform.org/alexa_seq/"><span class="underline">http://www.alexaplatform.org/alexa_seq/</span></a>.</td>
<td>1</td>
</tr>
<tr class="odd">
<td colspan="8">
<p><strong>h. Cell deconvolution</strong></p>
</td>

</tr>
<tr class="even">
<td>TIMER2.0<a href="https://paperpile.com/c/VjSfVH/J5Nw"><sup>134</sup></a></td>
<td>2020</td>
<td>Web server for comprehensive analysis of Tumor-Infiltrating Immune Cells.</td>
<td><p>Web-tool</p>
<p>R, Javascript</p></td>
<td>N/A</td>
<td>+</td>
<td><a href="https://github.com/taiwenli/TIMER"><span class="underline">https://github.com/taiwenli/TIMER</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>CIBERSORTx<a href="https://paperpile.com/c/VjSfVH/lg5HW"><sup>135</sup></a></td>
<td>2019</td>
<td><p>Impute gene expression profiles and provide an</p>
<p>estimation of the abundances of member cell types in a mixed cell population.</p></td>
<td><p>Web-tool</p>
<p>Java, R</p></td>
<td>N/A</td>
<td>+</td>
<td><a href="https://cibersortx.stanford.edu/"><span class="underline">https://cibersortx.stanford.edu/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>quanTIseq<a href="https://paperpile.com/c/VjSfVH/WG5s9"><sup>136</sup></a></td>
<td>2019</td>
<td>Quantify the fractions of ten immune cell types from bulk RNA-sequencing data.</td>
<td>R, Shell (Bash)</td>
<td>DockerHub</td>
<td>+</td>
<td><a href="https://urldefense.com/v3/__https://icbi.i-med.ac.at/quantiseq__;!!LIr3w8kk_Xxm!8I3jJ1Ohwgxxj7T2WzCPrAPGUzWpL74UMgqHVuJyzN3d_qimYh5begAcOYVGWVA$"><span class="underline">https://icbi.i-med.ac.at/quantiseq</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>Immunedeconv<a href="https://paperpile.com/c/VjSfVH/Q8hYH"><sup>137</sup></a></td>
<td>2019</td>
<td>Benchmarking of transcriptome-based cell-type quantification methods for immuno-oncology</td>
<td>R</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/icbi-lab/immunedeconv"><span class="underline">https://github.com/icbi-lab/immunedeconv</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Linseed<a href="https://paperpile.com/c/VjSfVH/OYixv"><sup>138</sup></a></td>
<td>2019</td>
<td>Deconvolution of cellular mixtures based on linearity of transcriptional signatures.</td>
<td>C++, R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/ctlab/LinSeed"><span class="underline">https://github.com/ctlab/LinSeed</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>deconvSEQ<a href="https://paperpile.com/c/VjSfVH/cCAZ9"><sup>139</sup></a></td>
<td>2019</td>
<td>Deconvolution of cell mixture distribution based on a generalized linear model.</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/rosedu1/deconvSeq">https://github.com/rosedu1/deconvSeq</a></td>
<td>1</td>
</tr>
<tr class="even">
<td>CDSeq<a href="https://paperpile.com/c/VjSfVH/bX8mF"><sup>140</sup></a></td>
<td>2019</td>
<td>Simultaneously estimate both cell-type proportions and cell-type-specific expression profiles.</td>
<td>MATLAB, R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/kkang7/CDSeq_R_Package"><span class="underline">https://github.com/kkang7/CDSeq_R_Package</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Dtangle<a href="https://paperpile.com/c/VjSfVH/V8Vl2"><sup>141</sup></a></td>
<td>2019</td>
<td>Estimates cell type proportions using publicly available, often cross-platform, reference data.</td>
<td>R</td>
<td>Anaconda,CRAN</td>
<td>++</td>
<td><a href="https://github.com/gjhunt/dtangle"><span class="underline">https://github.com/gjhunt/dtangle</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>GEDIT<a href="https://paperpile.com/c/VjSfVH/Q9P9g"><sup>142</sup></a></td>
<td>2019</td>
<td>Estimate cell type abundances.</td>
<td><p>Web based tool</p>
<p>Python, R</p></td>
<td>N/A</td>
<td>+</td>
<td><a href="http://webtools.mcdb.ucla.edu/"><span class="underline">http://webtools.mcdb.ucla.edu/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>SaVant<a href="https://paperpile.com/c/VjSfVH/svCH0"><sup>143</sup></a></td>
<td>2017</td>
<td>Web based tool for sample level visualization of molecular signatures in gene expression profiles.</td>
<td>Javascript, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://newpathways.mcdb.ucla.edu/savant"><span class="underline">http://newpathways.mcdb.ucla.edu/savant</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>EPIC<a href="https://paperpile.com/c/VjSfVH/EcMJl"><sup>144</sup></a></td>
<td>2017</td>
<td>Simultaneously estimates the fraction of cancer and immune cell types.</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://urldefense.com/v3/__http://epic.gfellerlab.org/__;!!LIr3w8kk_Xxm!_Ae4UYD0IkB4umt_LaL32HadyyK7IePPXv9wmbdb3WkVJtyBphUoMjbY-2dNwA4$"><span class="underline">http://epic.gfellerlab.org/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>WSCUnmix<a href="https://paperpile.com/c/VjSfVH/x3lNr"><sup>145</sup></a></td>
<td>2017</td>
<td>Automated deconvolution of structured mixtures.</td>
<td>MATLAB</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/tedroman/WSCUnmix"><span class="underline">https://github.com/tedroman/WSCUnmix</span></a></td>
<td><span class="underline">1</span></td>
</tr>
<tr class="even">
<td>Infino<a href="https://paperpile.com/c/VjSfVH/FpiwB"><sup>146</sup></a></td>
<td>2017</td>
<td>Deconvolves bulk RNA-seq into cell type abundances and captures gene expression variability in a Bayesian model to measure deconvolution uncertainty.</td>
<td>R, Python</td>
<td>Docker Hub</td>
<td>++</td>
<td><a href="https://github.com/hammerlab/infino"><span class="underline">https://github.com/hammerlab/infino</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>MCP-counter<a href="https://paperpile.com/c/VjSfVH/DcDBh"><sup>147</sup></a></td>
<td>2016</td>
<td>Estimating the population abundance of tissue-infiltrating immune and stromal cell populations.</td>
<td>R</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://github.com/ebecht/MCPcounter"><span class="underline">https://github.com/ebecht/MCPcounter</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>CellCode<a href="https://paperpile.com/c/VjSfVH/cHED9"><sup>148</sup></a></td>
<td>2015</td>
<td>Latent variable approach to differential expression analysis for heterogeneous cell populations.</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://www.pitt.edu/~mchikina/CellCODE/"><span class="underline">http://www.pitt.edu/~mchikina/CellCODE/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>PERT<a href="https://paperpile.com/c/VjSfVH/VDTaO"><sup>149</sup></a></td>
<td>2012</td>
<td>Probabilistic expression deconvolution method.</td>
<td>MATLAB</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/gquon/PERT"><span class="underline">https://github.com/gquon/PERT</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td colspan="8">
<p><strong>i. Immune repertoire profiling</strong></p></td>

</tr>
<tr class="odd">
<td>ImReP<a href="https://paperpile.com/c/VjSfVH/jYQwl"><sup>150</sup></a></td>
<td>2018</td>
<td>Profiling immunoglobulin repertoires across multiple human tissues.</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/mandricigor/imrep/wiki"><span class="underline">https://github.com/mandricigor/imrep/wiki</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>TRUST (T cell)<a href="https://paperpile.com/c/VjSfVH/lSZQ2"><sup>151</sup></a></td>
<td>2016</td>
<td>Landscape of tumor-infiltrating T cell repertoire of human cancers.</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/liulab-dfci/TRUST4"><span class="underline">https://github.com/liulab-dfci/TRUST4</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>V’DJer <a href="https://paperpile.com/c/VjSfVH/2Qvvf"><sup>152</sup></a></td>
<td>2016</td>
<td>Assembly-based inference of B-cell receptor repertoires from short reads with V’DJer.</td>
<td>C, C++</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/mozack/vdjer"><span class="underline">https://github.com/mozack/vdjer</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>IgBlast-based pipeline<a href="https://paperpile.com/c/VjSfVH/kps03"><sup>153</sup></a></td>
<td>2016</td>
<td>Statistical inference of a convergent antibody repertoire response.</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://www.ncbi.nlm.nih.gov/igblast/"><span class="underline">https://www.ncbi.nlm.nih.gov/igblast/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>MiXCR<a href="https://paperpile.com/c/VjSfVH/OBNcB"><sup>154</sup></a></td>
<td>2015</td>
<td>Processes big immunome data from raw sequences to quantitated clonotypes.</td>
<td>Java</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/milaboratory/mixcr"><span class="underline">https://github.com/milaboratory/mixcr</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td colspan="8">
<p><strong>j. Allele specific expression</strong></p></td>

</tr>
<tr class="odd">
<td>EAGLE<a href="https://paperpile.com/c/VjSfVH/k2fyX"><sup>155</sup></a></td>
<td>2017</td>
<td>Bayesian model for identifying GxE interactions based on associations between environmental variables and allele-specific expression.</td>
<td>C++, R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/davidaknowles/eagle"><span class="underline">https://github.com/davidaknowles/eagle</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>ANEVA-DOT/ANEVA<a href="https://paperpile.com/c/VjSfVH/9dQRu"><sup>156</sup></a></td>
<td>2019</td>
<td>Identify ASE outlier genes / Quantify genetic variation in gene dosage from ASE data.</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="https://github.com/PejLab/ANEVA-DOT"><span class="underline">https://github.com/PejLab/ANEVA-DOT</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>aFC<a href="https://paperpile.com/c/VjSfVH/ToouV"><sup>157</sup></a></td>
<td>2017</td>
<td>Quantifying the regulatory effect size of cis-acting genetic variation</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/secastel/aFC"><span class="underline">https://github.com/secastel/aFC</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>phASER<a href="https://paperpile.com/c/VjSfVH/EP329"><sup>158</sup></a></td>
<td>2016</td>
<td>Uses readback phasing to produce haplotype level ASE data (as opposed to SNP level)</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/secastel/phaser"><span class="underline">https://github.com/secastel/phaser</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>RASQUAL<a href="https://paperpile.com/c/VjSfVH/Nlke8"><sup>159</sup></a></td>
<td>2016</td>
<td>Maps QTLs for sequenced based cellular traits by combining population and allele-specific signals.</td>
<td>C, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/natsuhiko/rasqual"><span class="underline">https://github.com/natsuhiko/rasqual</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>allelecounter<a href="https://paperpile.com/c/VjSfVH/0eCk9"><sup>160</sup></a></td>
<td>2015</td>
<td>Generate ASE data from RNAseq data and a genotype file.</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/secastel/allelecounter"><span class="underline">https://github.com/secastel/allelecounter</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>WASP<a href="https://paperpile.com/c/VjSfVH/BF8iH"><sup>161</sup></a></td>
<td>2015</td>
<td>Unbiased allele-specific read mapping and discovery of molecular QTLs</td>
<td><p>C</p>
<p>Python</p></td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/bmvdgeijn/WASP/"><span class="underline">https://github.com/bmvdgeijn/WASP/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Mamba<a href="https://paperpile.com/c/VjSfVH/EjxUV"><sup>162</sup></a></td>
<td>2015</td>
<td>Compares different patterns of ASE across tissues</td>
<td>R</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://www.well.ox.ac.uk/~rivas/mamba/"><span class="underline">http://www.well.ox.ac.uk/~rivas/mamba/</span></a>.</td>
<td>2</td>
</tr>
<tr class="odd">
<td>MBASED <a href="https://paperpile.com/c/VjSfVH/qFoNp"><sup>163</sup></a></td>
<td>2014</td>
<td>Allele-specific expression detection in cancer tissues and cell lines</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/MBASED.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/MBASED.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Allim<a href="https://paperpile.com/c/VjSfVH/iWvro"><sup>164</sup></a></td>
<td>2013</td>
<td>Estimates allele-specific gene expression.</td>
<td>Python, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/allim/"><span class="underline">https://sourceforge.net/projects/allim/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>AlleleSeq<a href="https://paperpile.com/c/VjSfVH/gQ5lH"><sup>165</sup></a></td>
<td>2011</td>
<td>Identifies allele-specific events in mapped reads between maternal and paternal alleles.</td>
<td>Python, Shell</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://alleleseq.gersteinlab.org/"><span class="underline">http://alleleseq.gersteinlab.org/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td colspan="8">
<p><strong>k. Viral detection</strong></p></td>

</tr>
<tr class="odd">
<td>ROP<a href="https://paperpile.com/c/VjSfVH/9EvyR"><sup>166</sup></a></td>
<td>2018</td>
<td>Dumpster diving in RNA-sequencing to find the source of 1 trillion reads across diverse adult human tissues</td>
<td>Python, Shell (Bash)</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/smangul1/rop"><span class="underline">https://github.com/smangul1/rop</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>RNA CoMPASS<a href="https://paperpile.com/c/VjSfVH/rgUDB"><sup>167</sup></a></td>
<td>2014</td>
<td>Simultaneous analysis of transcriptomes and metatranscriptomes from diverse biological specimens.</td>
<td>Perl, Shell, Java</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://rnacompass.sourceforge.net/"><span class="underline">http://rnacompass.sourceforge.net/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>VirusSeq<a href="https://paperpile.com/c/VjSfVH/B4LYt"><sup>168</sup></a></td>
<td>2013</td>
<td>Identify viruses and their integration sites using next-generation sequencing of human cancer tissues</td>
<td>Perl, Shell (Bash)</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://odin.mdacc.tmc.edu/~xsu1/VirusSeq.html"><span class="underline">http://odin.mdacc.tmc.edu/∼xsu1/VirusSeq.htm</span></a>l</td>
<td>2</td>
</tr>
<tr class="even">
<td>VirusFinder<a href="https://paperpile.com/c/VjSfVH/OaVfr"><sup>169</sup></a></td>
<td>2013</td>
<td>Detection of Viruses and Their Integration Sites in Host Genomes through Next Generation Sequencing Data</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://bioinfo.mc.vanderbilt.edu/VirusFinder/"><span class="underline">http://bioinfo.mc.vanderbilt.edu/VirusFinder/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td colspan="8">
<p><strong>l. Fusion detection</strong></p>
</td>

</tr>
<tr class="even">
<td>INTEGRATE-Vis<a href="https://paperpile.com/c/VjSfVH/No9n"><sup>170</sup></a></td>
<td>2017</td>
<td><p>Generates plots focused on annotating each gene fusion at the transcript-</p>
<p>and protein-level and assessing expression across samples.</p></td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/ChrisMaherLab/INTEGRATE-Vis"><span class="underline">https://github.com/ChrisMaherLab/INTEGRATE-Vis</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>INTEGRATE-Neo<a href="https://paperpile.com/c/VjSfVH/Dbnw"><sup>171</sup></a></td>
<td>2017</td>
<td><p>Gene fusion neoantigen discovery tool, which uses RNA-Seq reads and is capable of</p>
<p>reporting tumor-specific peptides recognizable by immune molecules.</p></td>
<td>Python, C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/ChrisMaherLab/INTEGRATE-Neo"><span class="underline">https://github.com/ChrisMaherLab/INTEGRATE-Neo</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>INTEGRATE<a href="https://paperpile.com/c/VjSfVH/Ed4v"><sup>172</sup></a></td>
<td>2016</td>
<td>Capable of integrating aligned RNA-seq and WGS reads and characterizes the quality of predictions.</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/integrate-fusion/"><span class="underline">https://sourceforge.net/projects/integrate-fusion/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>TRUP<a href="https://paperpile.com/c/VjSfVH/ruXRb"><sup>173</sup></a></td>
<td>2015</td>
<td>Combines split-read and read-pair analysis with de novo regional assembly for the identification of chimeric transcripts in cancer specimens.</td>
<td>C++, Perl, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/ruping/TRUP"><span class="underline">https://github.com/ruping/TRUP</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>PRADA<a href="https://paperpile.com/c/VjSfVH/lWjhE"><sup>174</sup></a></td>
<td>2014</td>
<td>Detect gene fusions but also performs alignments, transcriptome quantification; mainly integrated genome/transcriptome read mapping.</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://sourceforge.net/projects/prada/"><span class="underline">http://sourceforge.net/projects/prada/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Pegasus<a href="https://paperpile.com/c/VjSfVH/RvAnG"><sup>175</sup></a></td>
<td>2014</td>
<td>Annotation and prediction of biologically functional gene fusion candidates.</td>
<td>Java, Perl, Python, Shell (Bash)</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/RabadanLab/Pegasus"><span class="underline">https://github.com/RabadanLab/Pegasus</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>FusionCatcher<a href="https://paperpile.com/c/VjSfVH/en5rB"><sup>176</sup></a></td>
<td>2014</td>
<td>Finding somatic fusion genes</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://urldefense.com/v3/__https://sourceforge.net/projects/fusioncatcher/__;!!LIr3w8kk_Xxm!9-fm4_qW1l6D5_Hnd3_hkmiNoHt9k3JhwjLcVJrdpp5hLZ6GIKngjsojBEDYFj0$"><span class="underline">https://sourceforge.net/projects/fusioncatcher/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>FusionQ<a href="https://paperpile.com/c/VjSfVH/QdIdG"><sup>177</sup></a></td>
<td>2013</td>
<td>Gene fusion detection and quantification from paired-end RNA-seq</td>
<td>C++, Perl, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.wakehealth.edu/CTSB/Software/Software.htm">http://www.wakehealth.edu/CTSB/Software/Software.htm</a></td>
<td>2</td>
</tr>
<tr class="even">
<td>Barnacle<a href="https://paperpile.com/c/VjSfVH/FumZk"><sup>178</sup></a></td>
<td>2013</td>
<td>Detecting and characterizing tandem duplications and fusions in de novo transcriptome assemblies</td>
<td>Python, Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.bcgsc.ca/platform/bioinfo/software/barnacle"><span class="underline">http://www.bcgsc.ca/platform/bioinfo/software/barnacle</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>Dissect<a href="https://paperpile.com/c/VjSfVH/bG5xT"><sup>179</sup></a></td>
<td>2012</td>
<td>Detection and characterization of structural alterations in transcribed sequences</td>
<td>C</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://dissect-trans.sourceforge.net/"><span class="underline">http://dissect-trans.sourceforge.net</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>BreakFusion<a href="https://paperpile.com/c/VjSfVH/ayKUC"><sup>180</sup></a></td>
<td>2012</td>
<td>Targeted assembly-based identification of gene fusions</td>
<td>C++, Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://bioinformatics.mdanderson.org/public-software/breakfusion/"><span class="underline">https://bioinformatics.mdanderson.org/public-software/breakfusion/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>EricScript<a href="https://paperpile.com/c/VjSfVH/VEFIX"><sup>181</sup></a></td>
<td>2012</td>
<td>Identification of gene fusion products in paired-end RNA-seq data.</td>
<td>Perl, R, Shell (Bash)</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://ericscript.sourceforge.net/"><span class="underline">http://ericscript.sourceforge.net</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Bellerophontes<a href="https://paperpile.com/c/VjSfVH/LKdvn"><sup>182</sup></a></td>
<td>2012</td>
<td>Chimeric transcripts discovery based on fusion model.</td>
<td>Java, Perl, Shell (Bash)</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://eda.polito.it/bellerophontes/"><span class="underline">http://eda.polito.it/bellerophontes/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>GFML<a href="https://paperpile.com/c/VjSfVH/x8Mju"><sup>183</sup></a></td>
<td>2012</td>
<td>Standard format for organizing and representing the significant features of gene fusion data.</td>
<td>XML</td>
<td>N/A</td>
<td>+</td>
<td><a href="http://code.google.com/p/gfml-prototype/"><span class="underline">http://code.google.com/p/gfml-prototype/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>FusionHunter<a href="https://paperpile.com/c/VjSfVH/1yU6B"><sup>184</sup></a></td>
<td>2011</td>
<td>Identifies fusion transcripts from transcriptional analysis.</td>
<td>C++</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://urldefense.com/v3/__https://github.com/ma-compbio/FusionHunter__;!!LIr3w8kk_Xxm!_9Tx482BhnmkWcVNOp6NGV6Nxr7pVGUUfdErDggLUumoBZ_PKQcdl7fGOtLCL1g$"><span class="underline">https://github.com/ma-compbio/FusionHunter</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>ChimeraScan<a href="https://paperpile.com/c/VjSfVH/XYjEV"><sup>185</sup></a></td>
<td>2011</td>
<td>Identifying chimeric transcription.</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://code.google.com/archive/p/chimerascan/downloads"><span class="underline">https://code.google.com/archive/p/chimerascan/downloads</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>TopHat-fusion<a href="https://paperpile.com/c/VjSfVH/mKjT3"><sup>186</sup></a></td>
<td>2011</td>
<td>Discovery of novel fusion transcripts.</td>
<td>C++, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://ccb.jhu.edu/software/tophat/fusion_index.shtml"><span class="underline">http://ccb.jhu.edu/software/tophat/fusion_index.shtml</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>deFuse<a href="https://paperpile.com/c/VjSfVH/G5ofQ"><sup>187</sup></a></td>
<td>2011</td>
<td>Fusion discovery in tumor RNA-seq data.</td>
<td>C++, Perl, R</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/amcpherson/defuse/blob/master/README.md"><span class="underline">https://github.com/amcpherson/defuse/blob/master/README.md</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td colspan="8">
<p><strong>m. Detecting circRNA</strong></p>
</td>

</tr>
<tr class="odd">
<td>CIRIquant<a href="https://paperpile.com/c/VjSfVH/WHE3"><sup>188</sup></a></td>
<td>2020</td>
<td>Accurate quantification and differential expression analysis of circRNAs.</td>
<td>Python</td>
<td>PyPI</td>
<td>++</td>
<td><a href="https://sourceforge.net/projects/ciri/files/CIRIquant/"><span class="underline">https://sourceforge.net/projects/ciri/files/CIRIquant</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>CIRI-vis<a href="https://paperpile.com/c/VjSfVH/hyGc"><sup>189</sup></a></td>
<td>2020</td>
<td>Visualization of circRNA structures.</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/ciri/files/CIRI-vis/"><span class="underline">https://sourceforge.net/projects/ciri/files/CIRI-vis</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>Ularcirc<a href="https://paperpile.com/c/VjSfVH/jD6k"><sup>190</sup></a></td>
<td>2019</td>
<td>Analysis and visualisation of canonical and back splice junctions.</td>
<td>R</td>
<td>Bioconductor</td>
<td>++</td>
<td><a href="https://github.com/VCCRI/Ularcirc"><span class="underline">https://github.com/VCCRI/Ularcirc</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>CLEAR<a href="https://paperpile.com/c/VjSfVH/Yuyn"><sup>191</sup></a></td>
<td>2019</td>
<td>Circular and Linear RNA expression analysis.</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/YangLab/CLEAR"><span class="underline">https://github.com/YangLab/CLEAR</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>CIRI-full<a href="https://paperpile.com/c/VjSfVH/1fZA"><sup>192</sup></a></td>
<td>2019</td>
<td>Reconstruct and quantify full-length circular RNAs.</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/ciri/files/CIRI-full/"><span class="underline">https://sourceforge.net/projects/ciri/files/CIRI-full</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>circAST<a href="https://paperpile.com/c/VjSfVH/q27T"><sup>193</sup></a></td>
<td>2019</td>
<td>Full-length assembly and quantification of alternatively spliced isoforms in Circular RNAs</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/xiaofengsong/CircAST"><span class="underline">https://github.com/xiaofengsong/CircAST</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>CIRI2<a href="https://paperpile.com/c/VjSfVH/hFQa"><sup>194</sup></a></td>
<td>2018</td>
<td>Denovo circRNA identification</td>
<td>Pearl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/ciri/files/CIRI2/"><span class="underline">https://sourceforge.net/projects/ciri/files/CIRI2</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Sailfish-cir<a href="https://paperpile.com/c/VjSfVH/drZh"><sup>195</sup></a></td>
<td>2017</td>
<td>Quantification of circRNAs using model-based framework</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/zerodel/sailfish-cir"><span class="underline">https://github.com/zerodel/sailfish-cir</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>CircComPara<a href="https://paperpile.com/c/VjSfVH/Yt5i"><sup>196</sup></a></td>
<td>2017</td>
<td>Multi‐method detection of circRNAs</td>
<td>R, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://github.com/egaffo/CirComPara"><span class="underline">http://github.com/egaffo/CirComPara</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>UROBORUS<a href="https://paperpile.com/c/VjSfVH/xlRs"><sup>197</sup></a></td>
<td>2016</td>
<td>Computationally identifying circRNAs from RNA-seq data</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/WGLab/UROBORUS/tree/master/bin"><span class="underline">https://github.com/WGLab/UROBORUS/tree/master/bin</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>PTESFinder<a href="https://paperpile.com/c/VjSfVH/TPEe"><sup>198</sup></a></td>
<td>2016</td>
<td>Identification of non-co-linear transcripts</td>
<td>Shell, Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/ptesfinder-v1/"><span class="underline">https://sourceforge.net/projects/ptesfinder-v1/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>NCLscan<a href="https://paperpile.com/c/VjSfVH/LPP0"><sup>199</sup></a></td>
<td>2016</td>
<td>identification of non-co-linear transcripts (fusion, trans-splicing and circular RNA)</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/TreesLab/NCLscan"><span class="underline">https://github.com/TreesLab/NCLscan</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>DCC<a href="https://paperpile.com/c/VjSfVH/O1Mu"><sup>200</sup></a></td>
<td>2016</td>
<td>Specific identification and quantification of circRNA</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/dieterich-lab/DCC"><span class="underline">https://github.com/dieterich-lab/DCC</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>CIRI-AS<a href="https://paperpile.com/c/VjSfVH/yHjv"><sup>201</sup></a></td>
<td>2016</td>
<td>Identification of internal structure and alternative splicing events in circRNA</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/ciri/files/CIRI-AS/"><span class="underline">https://sourceforge.net/projects/ciri/files/CIRI-AS</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>circTest<a href="https://paperpile.com/c/VjSfVH/O1Mu"><sup>200</sup></a></td>
<td>2016</td>
<td>Differential expression analysis and plotting of circRNAs</td>
<td>R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/dieterich-lab/CircTest"><span class="underline">https://github.com/dieterich-lab/CircTest</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>CIRCexplorer2<a href="https://paperpile.com/c/VjSfVH/FoV4"><sup>202</sup></a></td>
<td>2016</td>
<td>Annotation and de novo assembly of circRNAs</td>
<td>Python</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="https://github.com/YangLab/CIRCexplorer2"><span class="underline">https://github.com/YangLab/CIRCexplorer2</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>KNIFE <a href="https://paperpile.com/c/VjSfVH/3928"><sup>203</sup></a></td>
<td>2015</td>
<td>Statistically based detection of circular and linear isoforms from RNA-seq data</td>
<td>Perl, Python, R</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/lindaszabo/KNIFE"><span class="underline">https://github.com/lindaszabo/KNIFE</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>circRNA_finder<a href="https://paperpile.com/c/VjSfVH/uF4b"><sup>204</sup></a></td>
<td>2014</td>
<td>Identification of circRNAs from RNA-seq data</td>
<td>Perl</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://github.com/orzechoj/circRNA_finder"><span class="underline">https://github.com/orzechoj/circRNA_finder</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>find_circ<a href="https://paperpile.com/c/VjSfVH/Dw4f"><sup>205</sup></a></td>
<td>2013</td>
<td>Identification of circRNAs based on head-to-tail spliced sequencing reads</td>
<td>Python</td>
<td>PyPI</td>
<td>++</td>
<td><a href="https://github.com/marvin-jens/find_circ"><span class="underline">https://github.com/marvin-jens/find_circ</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td colspan="8">
<p><strong>n. Small RNA detection</strong></p>
</td>

</tr>
<tr class="odd">
<td>miRTrace<a href="https://paperpile.com/c/VjSfVH/NHns"><sup>206</sup></a></td>
<td>2018</td>
<td>Quality control of miRNA-seq data, identifies cross-species contamination.</td>
<td>Java</td>
<td>Anaconda</td>
<td>+++</td>
<td><a href="https://urldefense.com/v3/__https://github.com/friedlanderlab/mirtrace__;!!LIr3w8kk_Xxm!-Q6UWFNrpDD7q2N7AHJN5WAV_beS3YKz1KtPqXZDMKCQuOKr2Nt4bYMFdMxRHzs$"><span class="underline">https://github.com/friedlanderlab/mirtrace</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>sRNAbench<a href="https://paperpile.com/c/VjSfVH/LqiH"><sup>207</sup></a></td>
<td>2015</td>
<td>Expression profiling of small RNAs, prediction of novel microRNAs, analysis of isomiRs, genome mapping and read length statistics.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://bioinfo5.ugr.es/srnatoolbox/srnabench/"><span class="underline">https://bioinfo5.ugr.es/srnatoolbox/srnabench/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>sRNAde<a href="https://paperpile.com/c/VjSfVH/LqiH"><sup>207</sup></a></td>
<td>2015</td>
<td>Detection of differentially expressed small RNAs based on three programs.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://bioinfo5.ugr.es/srnatoolbox/srnade/"><span class="underline">https://bioinfo5.ugr.es/srnatoolbox/srnade/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>sRNAblast<a href="https://paperpile.com/c/VjSfVH/LqiH"><sup>207</sup></a></td>
<td>2015</td>
<td>Aimed to determine the origin of unmapped or unassigned reads by means of a blast search against several remote databases.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://bioinfo5.ugr.es/srnatoolbox/srnablast/"><span class="underline">https://bioinfo5.ugr.es/srnatoolbox/srnablast/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>miRNAconsTarget<a href="https://paperpile.com/c/VjSfVH/LqiH"><sup>207</sup></a></td>
<td>2015</td>
<td>Consensus target prediction on user provided input data.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://bioinfo5.ugr.es/srnatoolbox/amirconstarget/"><span class="underline">https://bioinfo5.ugr.es/srnatoolbox/amirconstarget/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>sRNAjBrowser<a href="https://paperpile.com/c/VjSfVH/LqiH"><sup>207</sup></a></td>
<td>2015</td>
<td>Visualization of sRNA expression data in a genome context.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://bioinfo5.ugr.es/srnatoolbox/srnajbrowser/"><span class="underline">https://bioinfo5.ugr.es/srnatoolbox/srnajbrowser/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>sRNAjBrowserDE<a href="https://paperpile.com/c/VjSfVH/LqiH"><sup>207</sup></a></td>
<td>2015</td>
<td>Visualization of differential expression as a function of read length in a genome context.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://bioinfo5.ugr.es/srnatoolbox/srnajbrowserde/"><span class="underline">https://bioinfo5.ugr.es/srnatoolbox/srnajbrowserde/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>ShortStack<a href="https://paperpile.com/c/VjSfVH/3POq"><sup>208</sup></a></td>
<td>2013</td>
<td>Analyzes reference-aligned sRNA-seq data and performs. comprehensive <em>de novo</em> annotation and quantification of the inferred sRNA genes.</td>
<td>Perl</td>
<td>Anaconda</td>
<td>+++</td>
<td><a href="https://github.com/MikeAxtell/ShortStack"><span class="underline">https://github.com/MikeAxtell/ShortStack</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>mirTools 2.0<a href="https://paperpile.com/c/VjSfVH/eVLg"><sup>209</sup></a></td>
<td>2013</td>
<td>Detect, identify and profile various types, functional annotation and differentially expressed sRNAs.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="http://www.wzgenomics.cn/mr2_dev/"><span class="underline">http://www.wzgenomics.cn/mr2_dev/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>UEA sRNA Workbench<a href="https://paperpile.com/c/VjSfVH/bglJ"><sup>210</sup></a></td>
<td>2012</td>
<td>Complete analysis of single or multiple-sample small RNA datasets.</td>
<td>Web based tool, C++, Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/srnaworkbench/"><span class="underline">https://sourceforge.net/projects/srnaworkbench/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>miRDeep2<a href="https://paperpile.com/c/VjSfVH/2oAc"><sup>211</sup></a></td>
<td>2011</td>
<td>Discovers known and novel miRNAs, quantifies miRNA expression.</td>
<td>Perl</td>
<td>Anaconda</td>
<td>+++</td>
<td><a href="https://urldefense.com/v3/__https://github.com/rajewsky-lab/mirdeep2__;!!LIr3w8kk_Xxm!-Q6UWFNrpDD7q2N7AHJN5WAV_beS3YKz1KtPqXZDMKCQuOKr2Nt4bYMFAhQ1p1Q$"><span class="underline">https://github.com/rajewsky-lab/mirdeep2</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>miRanalyzer<a href="https://paperpile.com/c/VjSfVH/Ub5I"><sup>212</sup></a></td>
<td>2011</td>
<td>Detection of known and prediction of new microRNAs in high-throughput sequencing experiments.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="http://bioinfo2.ugr.es/miRanalyzer/miRanalyzer.php"><span class="underline">http://bioinfo2.ugr.es/miRanalyzer/miRanalyzer.php</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>SeqBuster<a href="https://paperpile.com/c/VjSfVH/PJ7m"><sup>213</sup></a></td>
<td>2010</td>
<td>Provides an automatized pre-analysis for sequence annotation for analysing small RNA data from Illumina sequencing.</td>
<td>Web based tool</td>
<td>Anaconda</td>
<td>+</td>
<td><a href="http://estivill_lab.crg.es/seqbuster"><span class="underline">http://estivill_lab.crg.es/seqbuster</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>DARIO<a href="https://paperpile.com/c/VjSfVH/NEN6"><sup>214</sup></a></td>
<td>2010</td>
<td>Allows to study short read data and provides a wide range of analysis features, including quality control, read normalization, and quantification.</td>
<td>Web based tool</td>
<td>N/A</td>
<td>+</td>
<td><a href="http://dario.bioinf.uni-leipzig.de/index.py"><span class="underline">http://dario.bioinf.uni-leipzig.de/index.py</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td colspan="8">
<p><strong>o. Visualization tools</strong></p>
</td>

</tr>
<tr class="even">
<td>BEAVR<a href="https://paperpile.com/c/VjSfVH/lOig"><sup>215</sup></a></td>
<td>2020</td>
<td>Facilitates interactive analysis and exploration of RNA-seq data, allowing statistical testing and visualization of the table of differentially expressed genes obtained.</td>
<td>R</td>
<td>Docker Hub</td>
<td>++</td>
<td><a href="https://github.com/developerpiru/BEAVR"><span class="underline">https://github.com/developerpiru/BEAVR</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>coseq<a href="https://paperpile.com/c/VjSfVH/nURXp"><sup>216</sup></a></td>
<td>2018</td>
<td>Co-expression analysis of sequencing data</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://bioconductor.org/packages/release/bioc/html/coseq.html"><span class="underline">https://bioconductor.org/packages/release/bioc/html/coseq.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>ReadXplorer<a href="https://paperpile.com/c/VjSfVH/yDVCn"><sup>217</sup></a></td>
<td>2016</td>
<td>Read mapping analysis and visualization</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://www.uni-giessen.de/fbz/fb08/Inst/bioinformatik/software/ReadXplorer"><span class="underline">https://www.uni-giessen.de/fbz/fb08/Inst/bioinformatik/software/ReadXplorer</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>Integrated Genome Browser<a href="https://paperpile.com/c/VjSfVH/S3Os"><sup>218</sup></a></td>
<td>2016</td>
<td>An interactive tool for visually analyzing tiling array data and enables quantification of alternative splicing</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.bioviz.org/"><span class="underline">http://www.bioviz.org/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>Sashimi plots<a href="https://paperpile.com/c/VjSfVH/beuLd"><sup>219</sup></a></td>
<td>2015</td>
<td>Quantitative visualization comparison of exon usage</td>
<td>Python</td>
<td>N/A</td>
<td>++</td>
<td><a href="http://miso.readthedocs.org/en/fastmiso/sashimi.html"><span class="underline">http://miso.readthedocs.org/en/fastmiso/sashimi.html</span></a></td>
<td><span class="underline">2</span></td>
</tr>
<tr class="odd">
<td>ASTALAVISTA<a href="https://paperpile.com/c/VjSfVH/8PM74"><sup>220</sup></a></td>
<td>2015</td>
<td>Reports all alternative splicing events reflected by transcript annotations</td>
<td>Java</td>
<td>Anaconda</td>
<td>++</td>
<td><a href="http://astalavista.sammeth.net/"><span class="underline">http://astalavista.sammeth.net/</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>RNASeqBrowser<a href="https://paperpile.com/c/VjSfVH/ZTg03"><sup>221</sup></a></td>
<td>2015</td>
<td>Incorporates and extends the functionality of the UCSC genome browser</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://www.australianprostatecentre.org/research/software/rnaseqbrowser"><span class="underline">http://www.australianprostatecentre.org/research/software/rnaseqbrowser</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>SplicePlot<a href="https://paperpile.com/c/VjSfVH/F6BgU"><sup>222</sup></a></td>
<td>2014</td>
<td>Visualizing splicing quantitative trait loci</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://montgomerylab.stanford.edu/spliceplot/index.html"><span class="underline">http://montgomerylab.stanford.edu/spliceplot/index.html</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>RNASeqViewer<a href="https://paperpile.com/c/VjSfVH/9Ywn0"><sup>223</sup></a></td>
<td>2014</td>
<td>Compare gene expression and alternative splicing</td>
<td>Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="https://sourceforge.net/projects/rnaseqbrowser/"><span class="underline">https://sourceforge.net/projects/rnaseqbrowser/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>PrimerSeq<a href="https://paperpile.com/c/VjSfVH/0lcAx"><sup>224</sup></a></td>
<td>2014</td>
<td>Systematic design and visualization of RT-PCR primers using RNA seq data</td>
<td>Java, C++, Python</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://primerseq.sourceforge.net/"><span class="underline">http://primerseq.sourceforge.net/</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Epiviz<a href="https://paperpile.com/c/VjSfVH/PdrCV"><sup>225</sup></a></td>
<td>2014</td>
<td>Combining algorithmic-statistical analysis and interactive visualization</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="https://epiviz.github.io/"><span class="underline">https://epiviz.github.io/</span></a></td>
<td>1</td>
</tr>
<tr class="odd">
<td>RNAbrowse<a href="https://paperpile.com/c/VjSfVH/VoN5Q"><sup>226</sup></a></td>
<td>2014</td>
<td>RNA-seq De Novo Assembly Results Browser</td>
<td>N/A</td>
<td>N/A</td>
<td>N/A</td>
<td><a href="http://bioinfo.genotoul.fr/RNAbrowse"><span class="underline">http://bioinfo.genotoul.fr/RNAbrowse</span></a></td>
<td>2</td>
</tr>
<tr class="even">
<td>ZENBU<a href="https://paperpile.com/c/VjSfVH/KL5oc"><sup>227</sup></a></td>
<td>2014</td>
<td>Interactive visualization and analysis of large-scale sequencing datasets</td>
<td>C++, Javascript</td>
<td>N/A</td>
<td>+</td>
<td><a href="https://urldefense.com/v3/__https://fantom.gsc.riken.jp/zenbu/__;!!LIr3w8kk_Xxm!4ArDxJT0BHPuxZxS6xE6holcue6ZpaDFY8lTrjurWiZUh8AyEMnNeVdYfD_xOjs$"><span class="underline">https://fantom.gsc.riken.jp/zenbu/</span></a></td>
<td>2</td>
</tr>
<tr class="odd">
<td>CummeRbund<a href="https://paperpile.com/c/VjSfVH/bhdyI"><sup>53</sup></a></td>
<td>2012</td>
<td>Navigate through data produced from a Cuffdiff RNA-seq differential expression analysis</td>
<td>R</td>
<td>Anaconda, Bioconductor</td>
<td>++</td>
<td><a href="http://bioconductor.org/packages/devel/bioc/html/cummeRbund.html"><span class="underline">http://bioconductor.org/packages/devel/bioc/html/cummeRbund.html</span></a></td>
<td>1</td>
</tr>
<tr class="even">
<td>Splicing Viewer<a href="https://paperpile.com/c/VjSfVH/d5CT4"><sup>228</sup></a></td>
<td>2012</td>
<td>Visualization of splice junctions and alternative splicing</td>
<td>Java</td>
<td>N/A</td>
<td>+++</td>
<td><a href="http://bioinformatics.zj.cn/splicingviewer">http://bioinformatics.zj.cn/splicingviewer</a>.</td>
<td>2</td>
</tr>
</tbody>
</table>



**Table 1:** Landscape of current computational methods for RNA-seq
analysis. We categorized RNA-seq tools published from 2008 to 2020 based
on processes in the RNA-seq pipeline and workflow; starting with data
quality control, read alignment, gene annotations, transcriptome
assembly, transcriptome quantification, differential expression, RNA
splicing, cell deconvolution, immune repertoire profiling, allele
specific expression, viral detection, fusion detection, detecting
circRNA, small RNA detection, and visualization tools. The third column
("Notable Features") presents key functionalities and methods used. The
fourth column ("Programing Language") presents the interface mode (e.g.,
GUI, web-based, programming language). The fifth column ("Package
Manager") highlights if a package manager such as Anaconda,
Bioconductor, CRAN, Docker Hub, pip, or PyPI is available for the tool.
We designated the assumed expertise level with a +, ++, or +++ in the
sixth column ("Required Expertise"). A "+" represents little to no
required expertise which would be assigned to a GUI based/web interface
tool. "++" was assigned to tools that require R and/or multiple
programming languages and whose software is located on Anaconda,
Bioconductor, CRAN, Docker Hub, pip, or PyPI. "+++" was assigned to
tools that require expertise in languages such as C, C++, Java, Python,
Perl, or Shell (Bash) and may or may not have a package manager present.
For each tool, we provide the links where the published tool software
can be found and downloaded ("Software"). In the seventh column ("Type
of URL"), each tool was assigned a "1" for web services designed to host
source code or "2" for others (e.g personal and/or university web
services).

**Methods**
===========

###

### **Determine features of RNA-seq tools**

We compiled 235 RNA-seq tools published between 2008 and 2020 which have
varying purposes and capabilities based on the type of analysis one is
conducting or the biological questions one is answering. We have
considered 15 domains of RNA-seq analysis (data quality control, read
alignment, gene annotations, transcriptome assembly, transcriptome
quantification, differential expression, RNA splicing, cell
deconvolution, immune repertoire profiling, allele specific expression,
viral detection, fusion detection, detecting circRNA, small RNA
detection, and visualization tools). After assigning each tool a
category based on its area of application, we highlighted each tool's
("Notable Features"). These notable features encompassed a range of
functionalities: purpose of the tool, which features made the tool
unique or not unique within its category and the way in which the tool
would take RNA-seq data as an input and present an output. We documented
whether each tool was web-based or required one or many programming
languages for installation and/or utilization ("Programming Language").
In addition to the programming languages, we highlighted whether a
package manager (e.g., Anaconda, Bioconductor, CRAN, Docker Hub, pip, or
PyPI) was available. Based on the combination of which programming
language was required and which package manager was available for each
tool, we assessed the required expertise needed to be able to install or
run the tool. If the tool was a web-based tool with no package manager
available or a web-based tool along with programming languages and
package manager present, we assigned the tool the little to none
required expertise of "+". If the tool required only R as a programming
language or along with other programming languages and had a package
manager present, the tool was assigned a required expertise of "++". In
addition, if programming languages other than R were required, and a
package manager was present, the tool was also assigned a "++". Lastly,
for tools that required a single programming language (other than R) or
multiple programming languages required and presence of no package
manager were assigned a "+++" for the most required expertise. Each
published tool had a designated software link where the tool can be
downloaded and installed. Based on the type of URL of the software
links, we assigned the tools a "1" or "2". An assignment of "1" meant
that the tool's software was hosted on a web service designed to host
source code. An assignment of "2" meant that the tool's software was
hosted on other web services (e.g personal and/or university web
services).

**References**
==============



1. Kumar, G., Ertel, A., Feldman, G., Kupper, J. & Fortina, P. iSeqQC:
a tool for expression-based quality control in RNA sequencing. *BMC
Bioinformatics* **21**, 56 (2020). <div id="ref1"></div>


2. Hicks, S. C. *et al.* Smooth quantile normalization.
*Biostatistics* **19**, 185--198 (2018).<div id="ref2"></div>


3. Babraham Bioinformatics - FastQC A Quality Control tool for High
Throughput Sequence Data. <div id="ref3"></div>

4\. [Guo, Y. *et al.* Multi-perspective quality control of Illumina
exome sequencing data using QC3. *Genomics* **103**, 323--328
(2014).](http://paperpile.com/b/VjSfVH/oXXz9)

5\. [Anvar, S. Y. *et al.* Determining the quality and complexity of
next-generation sequencing data without a reference genome. *Genome
Biol.* **15**, 555 (2014).](http://paperpile.com/b/VjSfVH/Z32Ou)

6\. [Yang, X. *et al.* HTQC: a fast quality control toolkit for Illumina
sequencing data. *BMC Bioinformatics* **14**, 33
(2013).](http://paperpile.com/b/VjSfVH/ZdZOL)

7\. [Bolger, A. M., Lohse, M. & Usadel, B. Trimmomatic: a flexible
trimmer for Illumina sequence data. *Bioinformatics* **30**, 2114--2120
(2014).](http://paperpile.com/b/VjSfVH/euzDZ)

8\. [Jiang, H., Lei, R., Ding, S.-W. & Zhu, S. Skewer: a fast and
accurate adapter trimmer for next-generation sequencing paired-end
reads. *BMC Bioinformatics* **15**, 182
(2014).](http://paperpile.com/b/VjSfVH/AYNs8)

9\. [Dodt, M., Roehr, J., Ahmed, R. & Dieterich, C. FLEXBAR---Flexible
Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
*Biology* vol. 1 895--905 (2012).](http://paperpile.com/b/VjSfVH/ZxxWC)

10\. [Kroll, K. W. *et al.* Quality Control for RNA-Seq (QuaCRS): An
Integrated Quality Control Pipeline. *Cancer Inform.* **13**, 7--14
(2014).](http://paperpile.com/b/VjSfVH/9pl3i)

11\. [Cabanski, C. R. *et al.* BlackOPs: increasing confidence in
variant detection through mappability filtering. *Nucleic Acids Res.*
**41**, e178 (2013).](http://paperpile.com/b/VjSfVH/kl83M)

12\. [Wang, L., Wang, S. & Li, W. RSeQC: quality control of RNA-seq
experiments. *Bioinformatics* **28**, 2184--2185
(2012).](http://paperpile.com/b/VjSfVH/hF0OY)

13\. [DeLuca, D. S. *et al.* RNA-SeQC: RNA-seq metrics for quality
control and process optimization. *Bioinformatics* **28**, 1530--1532
(2012).](http://paperpile.com/b/VjSfVH/UXSG9)

14\. [Jones, D. C., Ruzzo, W. L., Peng, X. & Katze, M. G. A new approach
to bias correction in RNA-Seq. *Bioinformatics* **28**, 921--928
(2012).](http://paperpile.com/b/VjSfVH/x9hen)

15\. [Lassmann, T., Hayashizaki, Y. & Daub, C. O. SAMStat: monitoring
biases in next generation sequencing data. *Bioinformatics* **27**,
130--131 (2011).](http://paperpile.com/b/VjSfVH/iAYlH)

16\. [Li, H. *et al.* The Sequence Alignment/Map format and SAMtools.
*Bioinformatics* **25**, 2078--2079
(2009).](http://paperpile.com/b/VjSfVH/VYwJx)

17\. [Liu, B. *et al.* deSALT: fast and accurate long transcriptomic
read alignment with de Bruijn graph-based index.
doi:](http://paperpile.com/b/VjSfVH/W81bL)[10.1101/612176](http://dx.doi.org/10.1101/612176)[.](http://paperpile.com/b/VjSfVH/W81bL)

18\. [Boratyn, G. M., Thierry-Mieg, J., Thierry-Mieg, D., Busby, B. &
Madden, T. L. Magic-BLAST, an accurate DNA and RNA-seq aligner for long
and short reads.
doi:](http://paperpile.com/b/VjSfVH/rXBdj)[10.1101/390013](http://dx.doi.org/10.1101/390013)[.](http://paperpile.com/b/VjSfVH/rXBdj)

19\. [Li, H. Minimap2: pairwise alignment for nucleotide sequences.
*Bioinformatics* **34**, 3094--3100
(2018).](http://paperpile.com/b/VjSfVH/c2pYX)

20\. [Lin, H.-N. & Hsu, W.-L. DART: a fast and accurate RNA-seq mapper
with a partitioning strategy. *Bioinformatics* vol. 34 190--197
(2018).](http://paperpile.com/b/VjSfVH/ilhnw)

21\. [Kahles, A., Behr, J. & Rätsch, G. MMR: a tool for read
multi-mapper resolution. *Bioinformatics* **32**, 770--772
(2016).](http://paperpile.com/b/VjSfVH/3n9g)

22\. [Bonfert, T., Kirner, E., Csaba, G., Zimmer, R. & Friedel, C. C.
ContextMap 2: fast and accurate context-based RNA-seq mapping. *BMC
Bioinformatics* **16**, 122
(2015).](http://paperpile.com/b/VjSfVH/rrOR6)

23\. [Kim, D., Langmead, B. & Salzberg, S. L. HISAT: a fast spliced
aligner with low memory requirements. *Nat. Methods* **12**, 357--360
(2015).](http://paperpile.com/b/VjSfVH/X61GA)

24\. [Hoffmann, S. *et al.* A multi-split mapping algorithm for circular
RNA, splicing, trans-splicing and fusion detection. *Genome Biol.*
**15**, R34 (2014).](http://paperpile.com/b/VjSfVH/cBerW)

25\. [Butterfield, Y. S. *et al.* JAGuaR: junction alignments to genome
for RNA-seq reads. *PLoS One* **9**, e102398
(2014).](http://paperpile.com/b/VjSfVH/Nhgde)

26\. [Philippe, N., Salson, M., Commes, T. & Rivals, E. CRAC: an
integrated approach to the analysis of RNA-seq reads. *Genome Biol.*
**14**, R30 (2013).](http://paperpile.com/b/VjSfVH/oInun)

27\. [Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner.
*Bioinformatics* **29**, 15--21
(2013).](http://paperpile.com/b/VjSfVH/HcGrE)

28\. [Liao, Y., Smyth, G. K. & Shi, W. The Subread aligner: fast,
accurate and scalable read mapping by seed-and-vote. *Nucleic Acids
Res.* **41**, e108 (2013).](http://paperpile.com/b/VjSfVH/75Mhc)

29\. [Kim, D. *et al.* TopHat2: accurate alignment of transcriptomes in
the presence of insertions, deletions and gene fusions. *Genome Biol.*
**14**, R36 (2013).](http://paperpile.com/b/VjSfVH/v4QLu)

30\. [Hu, J., Ge, H., Newman, M. & Liu, K. OSA: a fast and accurate
alignment tool for RNA-Seq. *Bioinformatics* vol. 28 1933--1934
(2012).](http://paperpile.com/b/VjSfVH/6Rnd9)

31\. [Zhang, Y. *et al.* PASSion: a pattern growth algorithm-based
pipeline for splice junction detection in paired-end RNA-Seq data.
*Bioinformatics* **28**, 479--486
(2012).](http://paperpile.com/b/VjSfVH/sfG8k)

32\. [Grant, G. R. *et al.* Comparative analysis of RNA-Seq alignment
algorithms and the RNA-Seq unified mapper (RUM). *Bioinformatics*
**27**, 2518--2528 (2011).](http://paperpile.com/b/VjSfVH/FosW8)

33\. [Huang, S. *et al.* SOAPsplice: Genome-Wide ab initio Detection of
Splice Junctions from RNA-Seq Data. *Front. Genet.* **2**, 46
(2011).](http://paperpile.com/b/VjSfVH/cabdO)

34\. [Wang, K. *et al.* MapSplice: accurate mapping of RNA-seq reads for
splice junction discovery. *Nucleic Acids Res.* **38**, e178
(2010).](http://paperpile.com/b/VjSfVH/L2tcL)

35\. [Au, K. F., Jiang, H., Lin, L., Xing, Y. & Wong, W. H. Detection of
splice junctions from paired-end RNA-seq data by SpliceMap. *Nucleic
Acids Res.* **38**, 4570--4578
(2010).](http://paperpile.com/b/VjSfVH/yCVDK)

36\. [Bryant, D. W., Jr, Shen, R., Priest, H. D., Wong, W.-K. & Mockler,
T. C. Supersplat\--spliced RNA-seq alignment. *Bioinformatics* **26**,
1500--1505 (2010).](http://paperpile.com/b/VjSfVH/kY2xi)

37\. [Dimon, M. T., Sorber, K. & DeRisi, J. L. HMMSplicer: a tool for
efficient and sensitive discovery of known and novel splice junctions in
RNA-Seq data. *PLoS One* **5**, e13875
(2010).](http://paperpile.com/b/VjSfVH/JrZZK)

38\. [De Bona, F., Ossowski, S., Schneeberger, K. & Rätsch, G. Optimal
spliced alignments of short sequence reads. *Bioinformatics* **24**,
i174--80 (2008).](http://paperpile.com/b/VjSfVH/hGOfW)

39\. [Tardaguila, M. *et al.* SQANTI: extensive characterization of
long-read transcript sequences for quality control in full-length
transcriptome identification and quantification. *Genome Res.* (2018)
doi:](http://paperpile.com/b/VjSfVH/T0s0)[10.1101/gr.222976.117](http://dx.doi.org/10.1101/gr.222976.117)[.](http://paperpile.com/b/VjSfVH/T0s0)

40\. [Musacchia, F., Basu, S., Petrosino, G., Salvemini, M. & Sanges, R.
Annocript: a flexible pipeline for the annotation of transcriptomes able
to identify putative long noncoding RNAs. *Bioinformatics* **31**,
2199--2201 (2015).](http://paperpile.com/b/VjSfVH/sAdu3)

41\. [Gao, Y., Wang, J. & Zhao, F. CIRI: an efficient and unbiased
algorithm for de novo circular RNA identification. *Genome Biology* vol.
16 (2015).](http://paperpile.com/b/VjSfVH/jtYuD)

42\. [Amman, F. *et al.* TSSAR: TSS annotation regime for dRNA-seq data.
*BMC Bioinformatics* **15**, 89
(2014).](http://paperpile.com/b/VjSfVH/16V6S)

43\. [Tang, A. D. *et al.* Full-length transcript characterization of
SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of
retained introns.
doi:](http://paperpile.com/b/VjSfVH/TTIUT)[10.1101/410183](http://dx.doi.org/10.1101/410183)[.](http://paperpile.com/b/VjSfVH/TTIUT)

44\. [Shao, M. & Kingsford, C. Accurate assembly of transcripts through
phase-preserving graph decomposition. *Nat. Biotechnol.* **35**,
1167--1169 (2017).](http://paperpile.com/b/VjSfVH/vP2QA)

45\. [Song, L., Sabunciyan, S. & Florea, L. CLASS2: accurate and
efficient splice variant annotation from RNA-seq reads. *Nucleic Acids
Res.* **44**, e98 (2016).](http://paperpile.com/b/VjSfVH/l0gw)

46\. [Pertea, M. *et al.* StringTie enables improved reconstruction of a
transcriptome from RNA-seq reads. *Nat. Biotechnol.* **33**, 290--295
(2015).](http://paperpile.com/b/VjSfVH/hhaWZ)

47\. [Chang, Z. *et al.* Bridger: a new framework for de novo
transcriptome assembly using RNA-seq data. *Genome Biology* vol. 16
(2015).](http://paperpile.com/b/VjSfVH/dTtcj)

48\. [Maretty, L., Sibbesen, J. A. & Krogh, A. Bayesian transcriptome
assembly. *Genome Biol.* **15**, 501
(2014).](http://paperpile.com/b/VjSfVH/qvR6K)

49\. [Le, H.-S., Schulz, M. H., McCauley, B. M., Hinman, V. F. &
Bar-Joseph, Z. Probabilistic error correction for RNA sequencing.
*Nucleic Acids Res.* **41**, e109
(2013).](http://paperpile.com/b/VjSfVH/Ql98S)

50\. [Bao, E., Jiang, T. & Girke, T. BRANCH: boosting RNA-Seq assemblies
with partial or related genomic sequences. *Bioinformatics* vol. 29
1250--1259 (2013).](http://paperpile.com/b/VjSfVH/BY7aZ)

51\. [Chu, H.-T. *et al.* EBARDenovo: highly accurate de novo assembly
of RNA-Seq with efficient chimera-detection. *Bioinformatics* **29**,
1004--1010 (2013).](http://paperpile.com/b/VjSfVH/ZGu0v)

52\. [Schulz, M. H., Zerbino, D. R., Vingron, M. & Birney, E. Oases:
robust de novo RNA-seq assembly across the dynamic range of expression
levels. *Bioinformatics* **28**, 1086--1092
(2012).](http://paperpile.com/b/VjSfVH/ZrrJT)

53\. [Trapnell, C. *et al.* Differential gene and transcript expression
analysis of RNA-seq experiments with TopHat and Cufflinks. *Nat.
Protoc.* **7**, 562--578 (2012).](http://paperpile.com/b/VjSfVH/bhdyI)

54\. [Feng, J., Li, W. & Jiang, T. Inference of isoforms from short
sequence reads. *J. Comput. Biol.* **18**, 305--321
(2011).](http://paperpile.com/b/VjSfVH/S3Nqy)

55\. [Li, W., Feng, J. & Jiang, T. IsoLasso: a LASSO regression approach
to RNA-Seq based transcriptome assembly. *J. Comput. Biol.* **18**,
1693--1707 (2011).](http://paperpile.com/b/VjSfVH/0Y01n)

56\. [Grabherr, M. G. *et al.* Full-length transcriptome assembly from
RNA-Seq data without a reference genome. *Nat. Biotechnol.* **29**,
644--652 (2011).](http://paperpile.com/b/VjSfVH/Klb2k)

57\. [Robertson, G. *et al.* De novo assembly and analysis of RNA-seq
data. *Nat. Methods* **7**, 909--912
(2010).](http://paperpile.com/b/VjSfVH/Y2S4e)

58\. [Guttman, M. *et al.* Ab initio reconstruction of cell
type--specific transcriptomes in mouse reveals the conserved
multi-exonic structure of lincRNAs. *Nature Biotechnology* vol. 28
503--510 (2010).](http://paperpile.com/b/VjSfVH/NbxHf)

59\. [Wyman, D. *et al.* A technology-agnostic long-read analysis
pipeline for transcriptome discovery and quantification.
doi:](http://paperpile.com/b/VjSfVH/Y1mRt)[10.1101/672931](http://dx.doi.org/10.1101/672931)[.](http://paperpile.com/b/VjSfVH/Y1mRt)

60\. [Patro, R., Duggal, G., Love, M. I., Irizarry, R. A. & Kingsford,
C. Salmon provides fast and bias-aware quantification of transcript
expression. *Nat. Methods* **14**, 417--419
(2017).](http://paperpile.com/b/VjSfVH/xwlWH)

61\. [Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal
probabilistic RNA-seq quantification. *Nat. Biotechnol.* **34**,
525--527 (2016).](http://paperpile.com/b/VjSfVH/nUid2)

62\. [nanoporetech. nanoporetech/wub.
*GitHub*](http://paperpile.com/b/VjSfVH/ZSFmM)
<https://github.com/nanoporetech/wub>[.](http://paperpile.com/b/VjSfVH/ZSFmM)

63\. [Schmid, M. W. & Grossniklaus, U. Rcount: simple and flexible
RNA-Seq read counting. *Bioinformatics* **31**, 436--437
(2015).](http://paperpile.com/b/VjSfVH/gYga0)

64\. [Anders, S., Pyl, P. T. & Huber, W. HTSeq\--a Python framework to
work with high-throughput sequencing data. *Bioinformatics* **31**,
166--169 (2015).](http://paperpile.com/b/VjSfVH/dNzvY)

65\. [Lee, S., Seo, C. H., Alver, B. H., Lee, S. & Park, P. J. EMSAR:
estimation of transcript abundance from RNA-seq data by
mappability-based segmentation and reclustering. *BMC Bioinformatics*
**16**, 278 (2015).](http://paperpile.com/b/VjSfVH/PdqIQ)

66\. [Finotello, F. *et al.* Reducing bias in RNA sequencing data: a
novel approach to compute counts. *BMC Bioinformatics* vol. 15 S7
(2014).](http://paperpile.com/b/VjSfVH/kujz1)

67\. [Hashimoto, T. B., Edwards, M. D. & Gifford, D. K. Universal count
correction for high-throughput sequencing. *PLoS Comput. Biol.* **10**,
e1003494 (2014).](http://paperpile.com/b/VjSfVH/w2lKf)

68\. [Patro, R., Mount, S. M. & Kingsford, C. Sailfish enables
alignment-free isoform quantification from RNA-seq reads using
lightweight algorithms. *Nat. Biotechnol.* **32**, 462--464
(2014).](http://paperpile.com/b/VjSfVH/CfcCV)

69\. [Rossell, D., Stephan-Otto Attolini, C., Kroiss, M. & Stöcker, A.
QUANTIFYING ALTERNATIVE SPLICING FROM PAIRED-END RNA-SEQUENCING DATA.
*Ann. Appl. Stat.* **8**, 309--330
(2014).](http://paperpile.com/b/VjSfVH/8K1xC)

70\. [Mangul, S. *et al.* Transcriptome assembly and quantification from
Ion Torrent RNA-Seq data. *BMC Genomics* **15 Suppl 5**, S7
(2014).](http://paperpile.com/b/VjSfVH/rqlp)

71\. [Liao, Y., Smyth, G. K. & Shi, W. featureCounts: an efficient
general purpose program for assigning sequence reads to genomic
features. *Bioinformatics* **30**, 923--930
(2014).](http://paperpile.com/b/VjSfVH/7Qa39)

72\. [Behr, J. *et al.* MITIE: Simultaneous RNA-Seq-based transcript
identification and quantification in multiple samples. *Bioinformatics*
**29**, 2529--2538 (2013).](http://paperpile.com/b/VjSfVH/IBZ5)

73\. [Mezlini, A. M. *et al.* iReckon: simultaneous isoform discovery
and abundance estimation from RNA-seq data. *Genome Res.* **23**,
519--529 (2013).](http://paperpile.com/b/VjSfVH/e16fg)

74\. [Roberts, A. & Pachter, L. Streaming fragment assignment for
real-time analysis of sequencing experiments. *Nat. Methods* **10**,
71--73 (2013).](http://paperpile.com/b/VjSfVH/rKytO)

75\. [Glaus, P., Honkela, A. & Rattray, M. Identifying differentially
expressed transcripts from RNA-seq data with biological variation.
*Bioinformatics* **28**, 1721--1728
(2012).](http://paperpile.com/b/VjSfVH/q2362)

76\. [Du, J. *et al.* IQSeq: integrated isoform quantification analysis
based on next-generation sequencing. *PLoS One* **7**, e29175
(2012).](http://paperpile.com/b/VjSfVH/v6vlI)

77\. [Li, W. & Jiang, T. Transcriptome assembly and isoform expression
level estimation from biased RNA-Seq reads. *Bioinformatics* **28**,
2914--2921 (2012).](http://paperpile.com/b/VjSfVH/xC6tA)

78\. [Xu, G. *et al.* SAMMate: a GUI tool for processing short read
alignments in SAM/BAM format. *Source Code Biol. Med.* **6**, 2
(2011).](http://paperpile.com/b/VjSfVH/WDgEC)

79\. [Kim, H., Bi, Y., Pal, S., Gupta, R. & Davuluri, R. V. IsoformEx:
isoform level gene expression estimation using weighted non-negative
least squares from mRNA-Seq data. *BMC Bioinformatics* **12**, 305
(2011).](http://paperpile.com/b/VjSfVH/lQkln)

80\. [Nicolae, M., Mangul, S., Măndoiu, I. I. & Zelikovsky, A.
Estimation of alternative splicing isoform frequencies from RNA-Seq
data. *Algorithms Mol. Biol.* **6**, 9
(2011).](http://paperpile.com/b/VjSfVH/XJjxm)

81\. [Li, B. & Dewey, C. N. RSEM: accurate transcript quantification
from RNA-Seq data with or without a reference genome. *BMC
Bioinformatics* vol. 12 (2011).](http://paperpile.com/b/VjSfVH/YTV4V)

82\. [Risso, D., Schwartz, K., Sherlock, G. & Dudoit, S. GC-Content
Normalization for RNA-Seq Data. *BMC Bioinformatics* **12**, 1--17
(2011).](http://paperpile.com/b/VjSfVH/nZEoj)

83\. [Turro, E. *et al.* Haplotype and isoform specific expression
estimation using multi-mapping RNA-seq reads. *Genome Biol.* **12**, R13
(2011).](http://paperpile.com/b/VjSfVH/iRLTs)

84\. [Katz, Y., Wang, E. T., Airoldi, E. M. & Burge, C. B. Analysis and
design of RNA sequencing experiments for identifying isoform regulation.
*Nat. Methods* **7**, 1009--1015
(2010).](http://paperpile.com/b/VjSfVH/hDihL)

85\. [Richard, H. *et al.* Prediction of alternative isoforms from exon
expression levels in RNA-Seq experiments. *Nucleic Acids Res.* **38**,
e112 (2010).](http://paperpile.com/b/VjSfVH/GJ3dA)

86\. [Jiang, H. & Wong, W. H. Statistical inferences for isoform
expression in RNA-Seq. *Bioinformatics* **25**, 1026--1032
(2009).](http://paperpile.com/b/VjSfVH/nCaHC)

87\. [Bohnert, R., Behr, J. & Rätsch, G. Transcript quantification with
RNA-Seq data. *BMC Bioinformatics* vol. 10
(2009).](http://paperpile.com/b/VjSfVH/wuq2V)

88\. [Mortazavi, A., Williams, B. A., McCue, K., Schaeffer, L. & Wold,
B. Mapping and quantifying mammalian transcriptomes by RNA-Seq. *Nat.
Methods* **5**, 621--628 (2008).](http://paperpile.com/b/VjSfVH/O6Lej)

89\. [Zhu, A., Srivastava, A., Ibrahim, J. G., Patro, R. & Love, M. I.
Nonparametric expression analysis using inferential replicate counts.
*Nucleic Acids Res.* **47**, e105
(2019).](http://paperpile.com/b/VjSfVH/Tn9qR)

90\. [Gunady, M. K., Mount, S. M. & Corrada Bravo, H. Yanagi: Fast and
interpretable segment-based alternative splicing and gene expression
analysis. *BMC Bioinformatics* **20**, 421
(2019).](http://paperpile.com/b/VjSfVH/CAWU)

91\. [Sterne-Weiler, T., Weatheritt, R. J., Best, A. J., Ha, K. C. H. &
Blencowe, B. J. Efficient and Accurate Quantitative Profiling of
Alternative Splicing Patterns of Any Complexity on a Laptop. *Mol. Cell*
**72**, 187--200.e6 (2018).](http://paperpile.com/b/VjSfVH/s2eM)

92\. [Spurr, L. *et al.* ReQTL -- an allele-level measure of
variation-expression genomic relationships.
doi:](http://paperpile.com/b/VjSfVH/RWWwc)[10.1101/464206](http://dx.doi.org/10.1101/464206)[.](http://paperpile.com/b/VjSfVH/RWWwc)

93\. [Tapial, J. *et al.* An atlas of alternative splicing profiles and
functional associations reveals new regulatory programs and genes that
simultaneously express multiple major isoforms. *Genome Res.* **27**,
1759--1768 (2017).](http://paperpile.com/b/VjSfVH/mkVp)

94\. [Frazee, A. C. *et al.* Ballgown bridges the gap between
transcriptome assembly and expression analysis. *Nat. Biotechnol.*
**33**, 243--246 (2015).](http://paperpile.com/b/VjSfVH/fWFKP)

95\. [Law, C. W., Chen, Y., Shi, W. & Smyth, G. K. voom: Precision
weights unlock linear model analysis tools for RNA-seq read counts.
*Genome Biol.* **15**, R29 (2014).](http://paperpile.com/b/VjSfVH/BuZHd)

96\. [Shen, S. *et al.* rMATS: robust and flexible detection of
differential alternative splicing from replicate RNA-Seq data. *Proc.
Natl. Acad. Sci. U. S. A.* **111**, E5593--601
(2014).](http://paperpile.com/b/VjSfVH/0run)

97\. [Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold
change and dispersion for RNA-seq data with DESeq2. *Genome Biol.*
**15**, 550 (2014).](http://paperpile.com/b/VjSfVH/PKb3X)

98\. [Davidson, N. M. & Oshlack, A. Corset: enabling differential gene
expression analysis for de novoassembled transcriptomes. *Genome
Biology* vol. 15 (2014).](http://paperpile.com/b/VjSfVH/jAg8X)

99\. [Gu, J., Wang, X., Halakivi-Clarke, L., Clarke, R. & Xuan, J.
BADGE: a novel Bayesian model for accurate abundance quantification and
differential analysis of RNA-Seq data. *BMC Bioinformatics* **15 Suppl
9**, S6 (2014).](http://paperpile.com/b/VjSfVH/FU1Ma)

100\. [Soneson, C. compcodeR\--an R package for benchmarking
differential expression methods for RNA-seq data. *Bioinformatics* vol.
30 2517--2518 (2014).](http://paperpile.com/b/VjSfVH/AsPdk)

101\. [Rau, A., Marot, G. & Jaffrézic, F. Differential meta-analysis of
RNA-seq data from multiple studies. *BMC Bioinformatics* **15**, 91
(2014).](http://paperpile.com/b/VjSfVH/qlyC3)

102\. [Clark, N. R. *et al.* The characteristic direction: a geometrical
approach to identify differentially expressed genes. *BMC
Bioinformatics* **15**, 79 (2014).](http://paperpile.com/b/VjSfVH/DUkqa)

103\. [Rau, A., Gallopin, M., Celeux, G. & Jaffrézic, F. Data-based
filtering for replicated high-throughput transcriptome sequencing
experiments. *Bioinformatics* **29**, 2146--2152
(2013).](http://paperpile.com/b/VjSfVH/KkVO6)

104\. [Bi, Y. & Davuluri, R. V. NPEBseq: nonparametric empirical
bayesian-based procedure for differential expression analysis of RNA-seq
data. *BMC Bioinformatics* **14**, 262
(2013).](http://paperpile.com/b/VjSfVH/tePCt)

105\. [Leng, N. *et al.* EBSeq: an empirical Bayes hierarchical model
for inference in RNA-seq experiments. *Bioinformatics* **29**,
1035--1043 (2013).](http://paperpile.com/b/VjSfVH/DFbQ5)

106\. [Yu, D., Huber, W. & Vitek, O. Shrinkage estimation of dispersion
in Negative Binomial models for RNA-seq experiments with small sample
size. *Bioinformatics* **29**, 1275--1282
(2013).](http://paperpile.com/b/VjSfVH/2YaM9)

107\. [Trapnell, C. *et al.* Differential analysis of gene regulation at
transcript resolution with RNA-seq. *Nat. Biotechnol.* **31**, 46--53
(2013).](http://paperpile.com/b/VjSfVH/orKVr)

108\. [Li, J. & Tibshirani, R. Finding consistent patterns: a
nonparametric approach for identifying differential expression in
RNA-Seq data. *Stat. Methods Med. Res.* **22**, 519--536
(2013).](http://paperpile.com/b/VjSfVH/q2TVc)

109\. [Wang, W., Qin, Z., Feng, Z., Wang, X. & Zhang, X. Identifying
differentially spliced genes from two groups of RNA-seq samples. *Gene*
**518**, 164--170 (2013).](http://paperpile.com/b/VjSfVH/3Ks2)

110\. [Tarazona, S., García-Alcalde, F., Dopazo, J., Ferrer, A. &
Conesa, A. Differential expression in RNA-seq: a matter of depth.
*Genome Res.* **21**, 2213--2223
(2011).](http://paperpile.com/b/VjSfVH/lZsbu)

111\. [Robinson, M. D., McCarthy, D. J. & Smyth, G. K. edgeR: a
Bioconductor package for differential expression analysis of digital
gene expression data. *Bioinformatics* **26**, 139--140
(2010).](http://paperpile.com/b/VjSfVH/pnOy0)

112\. [Wang, L., Feng, Z., Wang, X., Wang, X. & Zhang, X. DEGseq: an R
package for identifying differentially expressed genes from RNA-seq
data. *Bioinformatics* vol. 26 136--138
(2010).](http://paperpile.com/b/VjSfVH/Z0J7)

113\. [Li, Y. I. *et al.* Annotation-free quantification of RNA splicing
using LeafCutter. *Nat. Genet.* **50**, 151--158
(2018).](http://paperpile.com/b/VjSfVH/munu8)

114\. [Green, C. J., Gazzara, M. R. & Barash, Y. MAJIQ-SPEL: web-tool to
interrogate classical and complex splicing variations from RNA-Seq data.
*Bioinformatics* **34**, 300--302
(2018).](http://paperpile.com/b/VjSfVH/2Y6cD)

115\. [Vaquero-Garcia, J. *et al.* A new view of transcriptome
complexity and regulation through the lens of local splicing variations.
*Elife* **5**, e11752 (2016).](http://paperpile.com/b/VjSfVH/b9zV)

116\. [Kahles, A., Ong, C. S., Zhong, Y. & Rätsch, G. SplAdder:
identification, quantification and testing of alternative splicing
events from RNA-Seq data. *Bioinformatics* **32**, 1840--1847
(2016).](http://paperpile.com/b/VjSfVH/zvhJm)

117\. [Pulyakhina, I. *et al.* SplicePie: a novel analytical approach
for the detection of alternative, non-sequential and recursive splicing.
*Nucleic Acids Res.* **43**, 11068
(2015).](http://paperpile.com/b/VjSfVH/9WudK)

118\. [Alamancos, G. P., Pagès, A., Trincado, J. L., Bellora, N. &
Eyras, E. Leveraging transcript quantification for fast computation of
alternative splicing profiles. *RNA* **21**, 1521--1531
(2015).](http://paperpile.com/b/VjSfVH/5JQWU)

119\. [Mudvari, P. *et al.* SNPlice: variants that modulate Intron
retention from RNA-sequencing data. *Bioinformatics* **31**, 1191--1198
(2015).](http://paperpile.com/b/VjSfVH/5NuGs)

120\. [Niu, L., Huang, W., Umbach, D. M. & Li, L. IUTA: a tool for
effectively detecting differential isoform usage from RNA-Seq data. *BMC
Genomics* **15**, 862 (2014).](http://paperpile.com/b/VjSfVH/bt7r8)

121\. [Kimes, P. K. *et al.* SigFuge: single gene clustering of RNA-seq
reveals differential isoform usage among cancer samples. *Nucleic Acids
Res.* **42**, e113 (2014).](http://paperpile.com/b/VjSfVH/dc2Lc)

122\. [Gatto, A. *et al.* FineSplice, enhanced splice junction detection
and quantification: a novel pipeline based on the assessment of diverse
RNA-Seq alignment solutions. *Nucleic Acids Res.* **42**, e71
(2014).](http://paperpile.com/b/VjSfVH/itVqO)

123\. [Hu, Y. *et al.* PennSeq: accurate isoform-specific gene
expression quantification in RNA-Seq by modeling non-uniform read
distribution. *Nucleic Acids Res.* **42**, e20
(2014).](http://paperpile.com/b/VjSfVH/c1nBL)

124\. [Bernard, E., Jacob, L., Mairal, J. & Vert, J.-P. Efficient RNA
isoform identification and quantification from RNA-Seq data with network
flows. *Bioinformatics* **30**, 2447--2455
(2014).](http://paperpile.com/b/VjSfVH/xFcS0)

125\. [Ye, Z. *et al.* Computational analysis reveals a correlation of
exon-skipping events with splicing, transcription and epigenetic
factors. *Nucleic Acids Res.* **42**, 2856--2869
(2014).](http://paperpile.com/b/VjSfVH/bg83O)

126\. [Vitting-Seerup, K., Porse, B. T., Sandelin, A. & Waage, J. E.
spliceR: An R package for classification of alternative splicing and
prediction of coding potential from RNA-seq data.
doi:](http://paperpile.com/b/VjSfVH/ze0g0)[10.7287/peerj.preprints.80](http://dx.doi.org/10.7287/peerj.preprints.80)[.](http://paperpile.com/b/VjSfVH/ze0g0)

127\. [Park, J. W., Tokheim, C., Shen, S. & Xing, Y. Identifying
differential alternative splicing events from RNA sequencing data using
RNASeq-MATS. *Methods Mol. Biol.* **1038**, 171--179
(2013).](http://paperpile.com/b/VjSfVH/gNSyJ)

128\. [Aschoff, M. *et al.* SplicingCompass: differential splicing
detection using RNA-seq data. *Bioinformatics* **29**, 1141--1148
(2013).](http://paperpile.com/b/VjSfVH/iR3p2)

129\. [Hu, Y. *et al.* DiffSplice: the genome-wide detection of
differential splicing events with RNA-seq. *Nucleic Acids Res.* **41**,
e39 (2013).](http://paperpile.com/b/VjSfVH/36q9H)

130\. [Anders, S., Reyes, A. & Huber, W. Detecting differential usage of
exons from RNA-seq data. *Genome Res.* **22**, 2008--2017
(2012).](http://paperpile.com/b/VjSfVH/TbBMF)

131\. [Ryan, M. C., Cleland, J., Kim, R., Wong, W. C. & Weinstein, J. N.
SpliceSeq: a resource for analysis and visualization of RNA-Seq data on
alternative splicing and its functional impacts. *Bioinformatics*
**28**, 2385--2387 (2012).](http://paperpile.com/b/VjSfVH/WUPzf)

132\. [Brooks, A. N. *et al.* Conservation of an RNA regulatory map
between Drosophila and mammals. *Genome Res.* **21**, 193--202
(2011).](http://paperpile.com/b/VjSfVH/q9K4)

133\. [Griffith, M. *et al.* Alternative expression analysis by RNA
sequencing. *Nat. Methods* **7**, 843--847
(2010).](http://paperpile.com/b/VjSfVH/CvwSl)

134\. [Li, T. *et al.* TIMER2.0 for analysis of tumor-infiltrating
immune cells. *Nucleic Acids Res.* (2020)
doi:](http://paperpile.com/b/VjSfVH/J5Nw)[10.1093/nar/gkaa407](http://dx.doi.org/10.1093/nar/gkaa407)[.](http://paperpile.com/b/VjSfVH/J5Nw)

135\. [Newman, A. M. *et al.* Determining cell type abundance and
expression from bulk tissues with digital cytometry. *Nat. Biotechnol.*
**37**, 773--782 (2019).](http://paperpile.com/b/VjSfVH/lg5HW)

136\. [Finotello, F. *et al.* Molecular and pharmacological modulators
of the tumor immune contexture revealed by deconvolution of RNA-seq
data. *Genome Med.* **11**, 34
(2019).](http://paperpile.com/b/VjSfVH/WG5s9)

137\. [Sturm, G. *et al.* Comprehensive evaluation of
transcriptome-based cell-type quantification methods for
immuno-oncology. *Bioinformatics* **35**, i436--i445
(2019).](http://paperpile.com/b/VjSfVH/Q8hYH)

138\. [Zaitsev, K., Bambouskova, M., Swain, A. & Artyomov, M. N.
Complete deconvolution of cellular mixtures based on linearity of
transcriptional signatures. *Nat. Commun.* **10**, 2209
(2019).](http://paperpile.com/b/VjSfVH/OYixv)

139\. [Du, R., Carey, V. & Weiss, S. T. deconvSeq: deconvolution of cell
mixture distribution in sequencing data. *Bioinformatics* **35**,
5095--5102 (2019).](http://paperpile.com/b/VjSfVH/cCAZ9)

140\. [Kang, K. *et al.* CDSeq: A novel complete deconvolution method
for dissecting heterogeneous samples using gene expression data. *PLoS
Comput. Biol.* **15**, e1007510
(2019).](http://paperpile.com/b/VjSfVH/bX8mF)

141\. [Hunt, G. J., Freytag, S., Bahlo, M. & Gagnon-Bartsch, J. A.
dtangle: accurate and robust cell type deconvolution. *Bioinformatics*
**35**, 2093--2099 (2019).](http://paperpile.com/b/VjSfVH/V8Vl2)

142\. [Nadel, B. *et al.* The Gene Expression Deconvolution Interactive
Tool (GEDIT): Accurate Cell Type Quantification from Gene Expression
Data.
doi:](http://paperpile.com/b/VjSfVH/Q9P9g)[10.1101/728493](http://dx.doi.org/10.1101/728493)[.](http://paperpile.com/b/VjSfVH/Q9P9g)

143\. [Lopez, D. *et al.* SaVanT: a web-based tool for the sample-level
visualization of molecular signatures in gene expression profiles. *BMC
Genomics* vol. 18 (2017).](http://paperpile.com/b/VjSfVH/svCH0)

144\. [Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E. &
Gfeller, D. Simultaneous enumeration of cancer and immune cell types
from bulk tumor gene expression data. *Elife* **6**,
(2017).](http://paperpile.com/b/VjSfVH/EcMJl)

145\. [Roman, T., Xie, L. & Schwartz, R. Automated deconvolution of
structured mixtures from heterogeneous tumor genomic data. *PLoS Comput.
Biol.* **13**, e1005815 (2017).](http://paperpile.com/b/VjSfVH/x3lNr)

146\. [Zaslavsky, M., Novik, J. B., Chang, E. & Hammerbacher, J. Infino:
a Bayesian hierarchical model improves estimates of immune infiltration
into tumor microenvironment.
doi:](http://paperpile.com/b/VjSfVH/FpiwB)[10.1101/221671](http://dx.doi.org/10.1101/221671)[.](http://paperpile.com/b/VjSfVH/FpiwB)

147\. [Becht, E. *et al.* Estimating the population abundance of
tissue-infiltrating immune and stromal cell populations using gene
expression. *Genome Biol.* **17**, 218
(2016).](http://paperpile.com/b/VjSfVH/DcDBh)

148\. [Chikina, M., Zaslavsky, E. & Sealfon, S. C. CellCODE: a robust
latent variable approach to differential expression analysis for
heterogeneous cell populations. *Bioinformatics* **31**, 1584--1591
(2015).](http://paperpile.com/b/VjSfVH/cHED9)

149\. [Qiao, W. *et al.* PERT: A Method for Expression Deconvolution of
Human Blood Samples from Varied Microenvironmental and Developmental
Conditions. *PLoS Computational Biology* vol. 8 e1002838
(2012).](http://paperpile.com/b/VjSfVH/VDTaO)

150\. [Mangul, S. *et al.* Profiling immunoglobulin repertoires across
multiple human tissues by RNA Sequencing.
doi:](http://paperpile.com/b/VjSfVH/jYQwl)[10.1101/089235](http://dx.doi.org/10.1101/089235)[.](http://paperpile.com/b/VjSfVH/jYQwl)

151\. [Li, B. *et al.* Landscape of tumor-infiltrating T cell repertoire
of human cancers. *Nature Genetics* vol. 48 725--732
(2016).](http://paperpile.com/b/VjSfVH/lSZQ2)

152\. [Mose, L. E. *et al.* Assembly-based inference of B-cell receptor
repertoires from short read RNA sequencing data with V'DJer.
*Bioinformatics* vol. 32 3729--3734
(2016).](http://paperpile.com/b/VjSfVH/2Qvvf)

153\. [Strauli, N. B. & Hernandez, R. D. Statistical inference of a
convergent antibody repertoire response to influenza vaccine. *Genome
Med.* **8**, 60 (2016).](http://paperpile.com/b/VjSfVH/kps03)

154\. [Bolotin, D. A. *et al.* MiXCR: software for comprehensive
adaptive immunity profiling. *Nat. Methods* **12**, 380--381
(2015).](http://paperpile.com/b/VjSfVH/OBNcB)

155\. [Knowles, D. A. *et al.* Allele-specific expression reveals
interactions between genetic variation and environment. *Nat. Methods*
**14**, 699--702 (2017).](http://paperpile.com/b/VjSfVH/k2fyX)

156\. [Mohammadi, P. *et al.* Genetic regulatory variation in
populations informs transcriptome analysis in rare disease. *Science*
**366**, 351--356 (2019).](http://paperpile.com/b/VjSfVH/9dQRu)

157\. [Mohammadi, P., Castel, S. E., Brown, A. A. & Lappalainen, T.
Quantifying the regulatory effect size of cis-acting genetic variation
using allelic fold change.
doi:](http://paperpile.com/b/VjSfVH/ToouV)[10.1101/078717](http://dx.doi.org/10.1101/078717)[.](http://paperpile.com/b/VjSfVH/ToouV)

158\. [Castel, S. E., Mohammadi, P., Chung, W. K., Shen, Y. &
Lappalainen, T. Rare variant phasing and haplotypic expression from RNA
sequencing with phASER. *Nat. Commun.* **7**, 12817
(2016).](http://paperpile.com/b/VjSfVH/EP329)

159\. [Kumasaka, N., Knights, A. J. & Gaffney, D. J. Fine-mapping
cellular QTLs with RASQUAL and ATAC-seq. *Nat. Genet.* **48**, 206--213
(2016).](http://paperpile.com/b/VjSfVH/Nlke8)

160\. [Castel, S. E., Levy-Moonshine, A., Mohammadi, P., Banks, E. &
Lappalainen, T. Tools and best practices for data processing in allelic
expression analysis. *Genome Biol.* **16**, 195
(2015).](http://paperpile.com/b/VjSfVH/0eCk9)

161\. [van de Geijn, B., McVicker, G., Gilad, Y. & Pritchard, J. K.
WASP: allele-specific software for robust molecular quantitative trait
locus discovery. *Nat. Methods* **12**, 1061--1063
(2015).](http://paperpile.com/b/VjSfVH/BF8iH)

162\. [Pirinen, M. *et al.* Assessing allele-specific expression across
multiple tissues from RNA-seq read data. *Bioinformatics* **31**,
2497--2504 (2015).](http://paperpile.com/b/VjSfVH/EjxUV)

163\. [Mayba, O. *et al.* MBASED: allele-specific expression detection
in cancer tissues and cell lines. *Genome Biol.* **15**, 405
(2014).](http://paperpile.com/b/VjSfVH/qFoNp)

164\. [Pandey, R. V., Franssen, S. U., Futschik, A. & Schlötterer, C.
Allelic imbalance metre (Allim), a new tool for measuring
allele-specific gene expression with RNA-seq data. *Mol. Ecol. Resour.*
**13**, 740--745 (2013).](http://paperpile.com/b/VjSfVH/iWvro)

165\. [Rozowsky, J. *et al.* AlleleSeq: analysis of allele‐specific
expression and binding in a network framework. *Molecular Systems
Biology* vol. 7 522 (2011).](http://paperpile.com/b/VjSfVH/gQ5lH)

166\. [Mangul, S. *et al.* ROP: dumpster diving in RNA-sequencing to
find the source of 1 trillion reads across diverse adult human tissues.
*Genome Biol.* **19**, 36 (2018).](http://paperpile.com/b/VjSfVH/9EvyR)

167\. [Xu, G. *et al.* RNA CoMPASS: a dual approach for pathogen and
host transcriptome analysis of RNA-seq datasets. *PLoS One* **9**,
e89445 (2014).](http://paperpile.com/b/VjSfVH/rgUDB)

168\. [Chen, Y. *et al.* VirusSeq: software to identify viruses and
their integration sites using next-generation sequencing of human cancer
tissue. *Bioinformatics* **29**, 266--267
(2013).](http://paperpile.com/b/VjSfVH/B4LYt)

169\. [Wang, Q., Jia, P. & Zhao, Z. VirusFinder: software for efficient
and accurate detection of viruses and their integration sites in host
genomes through next generation sequencing data. *PLoS One* **8**,
e64465 (2013).](http://paperpile.com/b/VjSfVH/OaVfr)

170\. [Zhang, J., Gao, T. & Maher, C. A. INTEGRATE-Vis: a tool for
comprehensive gene fusion visualization. *Sci. Rep.* **7**, 17808
(2017).](http://paperpile.com/b/VjSfVH/No9n)

171\. [Zhang, J., Mardis, E. R. & Maher, C. A. INTEGRATE-neo: a pipeline
for personalized gene fusion neoantigen discovery. *Bioinformatics*
**33**, 555--557 (2017).](http://paperpile.com/b/VjSfVH/Dbnw)

172\. [Zhang, J. *et al.* INTEGRATE: gene fusion discovery using whole
genome and transcriptome data. *Genome Res.* **26**, 108--118
(2016).](http://paperpile.com/b/VjSfVH/Ed4v)

173\. [Fernandez-Cuesta, L. *et al.* Identification of novel fusion
genes in lung cancer using breakpoint assembly of transcriptome
sequencing data. *Genome Biol.* **16**, 7
(2015).](http://paperpile.com/b/VjSfVH/ruXRb)

174\. [Torres-García, W. *et al.* PRADA: pipeline for RNA sequencing
data analysis. *Bioinformatics* **30**, 2224--2226
(2014).](http://paperpile.com/b/VjSfVH/lWjhE)

175\. [Abate, F. *et al.* Pegasus: a comprehensive annotation and
prediction tool for detection of driver gene fusions in cancer. *BMC
Syst. Biol.* **8**, 97 (2014).](http://paperpile.com/b/VjSfVH/RvAnG)

176\. [Nicorici, D. *et al.* FusionCatcher - a tool for finding somatic
fusion genes in paired-end RNA-sequencing data.
doi:](http://paperpile.com/b/VjSfVH/en5rB)[10.1101/011650](http://dx.doi.org/10.1101/011650)[.](http://paperpile.com/b/VjSfVH/en5rB)

177\. [Liu, C., Ma, J., Chang, C. J. & Zhou, X. FusionQ: a novel
approach for gene fusion detection and quantification from paired-end
RNA-Seq. *BMC Bioinformatics* **14**, 193
(2013).](http://paperpile.com/b/VjSfVH/QdIdG)

178\. [Swanson, L. *et al.* Barnacle: detecting and characterizing
tandem duplications and fusions in transcriptome assemblies. *BMC
Genomics* **14**, 550 (2013).](http://paperpile.com/b/VjSfVH/FumZk)

179\. [Yorukoglu, D. *et al.* Dissect: detection and characterization of
novel structural alterations in transcribed sequences. *Bioinformatics*
**28**, i179--87 (2012).](http://paperpile.com/b/VjSfVH/bG5xT)

180\. [Chen, K. *et al.* BreakFusion: targeted assembly-based
identification of gene fusions in whole transcriptome paired-end
sequencing data. *Bioinformatics* **28**, 1923--1924
(2012).](http://paperpile.com/b/VjSfVH/ayKUC)

181\. [Benelli, M. *et al.* Discovering chimeric transcripts in
paired-end RNA-seq data by using EricScript. *Bioinformatics* **28**,
3232--3239 (2012).](http://paperpile.com/b/VjSfVH/VEFIX)

182\. [Abate, F. *et al.* Bellerophontes: an RNA-Seq data analysis
framework for chimeric transcripts discovery based on accurate fusion
model. *Bioinformatics* **28**, 2114--2121
(2012).](http://paperpile.com/b/VjSfVH/LKdvn)

183\. [Kalyana-Sundaram, S., Shanmugam, A. & Chinnaiyan, A. M. Gene
Fusion Markup Language: a prototype for exchanging gene fusion data.
*BMC Bioinformatics* vol. 13 269
(2012).](http://paperpile.com/b/VjSfVH/x8Mju)

184\. [Li, Y., Chien, J., Smith, D. I. & Ma, J. FusionHunter:
identifying fusion transcripts in cancer using paired-end RNA-seq.
*Bioinformatics* **27**, 1708--1710
(2011).](http://paperpile.com/b/VjSfVH/1yU6B)

185\. [Iyer, M. K., Chinnaiyan, A. M. & Maher, C. A. ChimeraScan: a tool
for identifying chimeric transcription in sequencing data.
*Bioinformatics* **27**, 2903--2904
(2011).](http://paperpile.com/b/VjSfVH/XYjEV)

186\. [Kim, D. & Salzberg, S. L. TopHat-Fusion: an algorithm for
discovery of novel fusion transcripts. *Genome Biol.* **12**, R72
(2011).](http://paperpile.com/b/VjSfVH/mKjT3)

187\. [McPherson, A. *et al.* deFuse: an algorithm for gene fusion
discovery in tumor RNA-Seq data. *PLoS Comput. Biol.* **7**, e1001138
(2011).](http://paperpile.com/b/VjSfVH/G5ofQ)

188\. [Zhang, J., Chen, S., Yang, J. & Zhao, F. Accurate quantification
of circular RNAs identifies extensive circular isoform switching events.
*Nat. Commun.* **11**, 90 (2020).](http://paperpile.com/b/VjSfVH/WHE3)

189\. [Zheng, Y. & Zhao, F. Visualization of circular RNAs and their
internal splicing events from transcriptomic data. *Bioinformatics*
**36**, 2934--2935 (2020).](http://paperpile.com/b/VjSfVH/hyGc)

190\. [Humphreys, D. T., Fossat, N., Demuth, M., Tam, P. P. L. & Ho, J.
W. K. Ularcirc: visualization and enhanced analysis of circular RNAs via
back and canonical forward splicing. *Nucleic Acids Res.* **47**, e123
(2019).](http://paperpile.com/b/VjSfVH/jD6k)

191\. [Ma, X.-K. *et al.* A CLEAR pipeline for direct comparison of
circular and linear RNA expression.
doi:](http://paperpile.com/b/VjSfVH/Yuyn)[10.1101/668657](http://dx.doi.org/10.1101/668657)[.](http://paperpile.com/b/VjSfVH/Yuyn)

192\. [Zheng, Y., Ji, P., Chen, S., Hou, L. & Zhao, F. Reconstruction of
full-length circular RNAs enables isoform-level quantification. *Genome
Med.* **11**, 2 (2019).](http://paperpile.com/b/VjSfVH/1fZA)

193\. [Wu, J. *et al.* CircAST: Full-length Assembly and Quantification
of Alternatively Spliced Isoforms in Circular RNAs. *Genomics Proteomics
Bioinformatics* **17**, 522--534
(2019).](http://paperpile.com/b/VjSfVH/q27T)

194\. [Gao, Y., Zhang, J. & Zhao, F. Circular RNA identification based
on multiple seed matching. *Brief. Bioinform.* **19**, 803--810
(2018).](http://paperpile.com/b/VjSfVH/hFQa)

195\. [Li, M. *et al.* Quantifying circular RNA expression from RNA-seq
data using model-based framework. *Bioinformatics* **33**, 2131--2139
(2017).](http://paperpile.com/b/VjSfVH/drZh)

196\. [Gaffo, E., Bonizzato, A., Kronnie, G. & Bortoluzzi, S.
CirComPara: A Multi‐Method Comparative Bioinformatics Pipeline to Detect
and Study circRNAs from RNA‐seq Data. *Non-Coding RNA* vol. 3 8
(2017).](http://paperpile.com/b/VjSfVH/Yt5i)

197\. [Song, X. *et al.* Circular RNA profile in gliomas revealed by
identification tool UROBORUS. *Nucleic Acids Res.* **44**, e87
(2016).](http://paperpile.com/b/VjSfVH/xlRs)

198\. [Izuogu, O. G. *et al.* PTESFinder: a computational method to
identify post-transcriptional exon shuffling (PTES) events. *BMC
Bioinformatics* **17**, 31 (2016).](http://paperpile.com/b/VjSfVH/TPEe)

199\. [Chuang, T.-J. *et al.* NCLscan: accurate identification of
non-co-linear transcripts (fusion, trans-splicing and circular RNA) with
a good balance between sensitivity and precision. *Nucleic Acids Res.*
**44**, e29 (2016).](http://paperpile.com/b/VjSfVH/LPP0)

200\. [Cheng, J., Metge, F. & Dieterich, C. Specific identification and
quantification of circular RNAs from sequencing data. *Bioinformatics*
**32**, 1094--1096 (2016).](http://paperpile.com/b/VjSfVH/O1Mu)

201\. [Gao, Y. *et al.* Comprehensive identification of internal
structure and alternative splicing events in circular RNAs. *Nat.
Commun.* **7**, 12060 (2016).](http://paperpile.com/b/VjSfVH/yHjv)

202\. [Zhang, X.-O. *et al.* Diverse alternative back-splicing and
alternative splicing landscape of circular RNAs. *Genome Res.* **26**,
1277--1287 (2016).](http://paperpile.com/b/VjSfVH/FoV4)

203\. [Szabo, L. *et al.* Statistically based splicing detection reveals
neural enrichment and tissue-specific induction of circular RNA during
human fetal development. *Genome Biol.* **16**, 126
(2015).](http://paperpile.com/b/VjSfVH/3928)

204\. [Westholm, J. O. *et al.* Genome-wide analysis of drosophila
circular RNAs reveals their structural and sequence properties and
age-dependent neural accumulation. *Cell Rep.* **9**, 1966--1980
(2014).](http://paperpile.com/b/VjSfVH/uF4b)

205\. [Memczak, S. *et al.* Circular RNAs are a large class of animal
RNAs with regulatory potency. *Nature* **495**, 333--338
(2013).](http://paperpile.com/b/VjSfVH/Dw4f)

206\. [Kang, W. *et al.* miRTrace reveals the organismal origins of
microRNA sequencing data. *Genome Biol.* **19**, 213
(2018).](http://paperpile.com/b/VjSfVH/NHns)

207\. [Rueda, A. *et al.* sRNAtoolbox: an integrated collection of small
RNA research tools. *Nucleic Acids Res.* **43**, W467--73
(2015).](http://paperpile.com/b/VjSfVH/LqiH)

208\. [Axtell, M. J. ShortStack: comprehensive annotation and
quantification of small RNA genes. *RNA* **19**, 740--751
(2013).](http://paperpile.com/b/VjSfVH/3POq)

209\. [Wu, J. *et al.* mirTools 2.0 for non-coding RNA discovery,
profiling, and functional annotation based on high-throughput
sequencing. *RNA Biol.* **10**, 1087--1092
(2013).](http://paperpile.com/b/VjSfVH/eVLg)

210\. [Stocks, M. B. *et al.* The UEA sRNA workbench: a suite of tools
for analysing and visualizing next generation sequencing microRNA and
small RNA datasets. *Bioinformatics* **28**, 2059--2061
(2012).](http://paperpile.com/b/VjSfVH/bglJ)

211\. [Friedländer, M. R., Mackowiak, S. D., Li, N., Chen, W. &
Rajewsky, N. miRDeep2 accurately identifies known and hundreds of novel
microRNA genes in seven animal clades. *Nucleic Acids Res.* **40**,
37--52 (2012).](http://paperpile.com/b/VjSfVH/2oAc)

212\. [Hackenberg, M., Rodríguez-Ezpeleta, N. & Aransay, A. M.
miRanalyzer: an update on the detection and analysis of microRNAs in
high-throughput sequencing experiments. *Nucleic Acids Res.* **39**,
W132--8 (2011).](http://paperpile.com/b/VjSfVH/Ub5I)

213\. [Pantano, L., Estivill, X. & Martí, E. SeqBuster, a bioinformatic
tool for the processing and analysis of small RNAs datasets, reveals
ubiquitous miRNA modifications in human embryonic cells. *Nucleic Acids
Res.* **38**, e34 (2010).](http://paperpile.com/b/VjSfVH/PJ7m)

214\. [Fasold, M., Langenberger, D., Binder, H., Stadler, P. F. &
Hoffmann, S. DARIO: a ncRNA detection and analysis tool for
next-generation sequencing experiments. *Nucleic Acids Research* vol. 39
W112--W117 (2011).](http://paperpile.com/b/VjSfVH/NEN6)

215\. [Perampalam, P. & Dick, F. A. BEAVR: a browser-based tool for the
exploration and visualization of RNA-seq data. *BMC Bioinformatics*
**21**, 221 (2020).](http://paperpile.com/b/VjSfVH/lOig)

216\. [Rau, A. & Maugis-Rabusseau, C. Transformation and model choice
for RNA-seq co-expression analysis. *Brief. Bioinform.* **19**, 425--436
(2018).](http://paperpile.com/b/VjSfVH/nURXp)

217\. [Hilker, R. *et al.* ReadXplorer 2---detailed read mapping
analysis and visualization from one single source. *Bioinformatics* vol.
32 3702--3708 (2016).](http://paperpile.com/b/VjSfVH/yDVCn)

218\. [Freese, N. H., Norris, D. C. & Loraine, A. E. Integrated genome
browser: visual analytics platform for genomics. *Bioinformatics*
**32**, 2089--2095 (2016).](http://paperpile.com/b/VjSfVH/S3Os)

219\. [Katz, Y. *et al.* Quantitative visualization of alternative exon
expression from RNA-seq data. *Bioinformatics* **31**, 2400--2402
(2015).](http://paperpile.com/b/VjSfVH/beuLd)

220\. [Foissac, S. & Sammeth, M. Analysis of alternative splicing events
in custom gene datasets by AStalavista. *Methods Mol. Biol.* **1269**,
379--392 (2015).](http://paperpile.com/b/VjSfVH/8PM74)

221\. [An, J. *et al.* RNASeqBrowser: a genome browser for simultaneous
visualization of raw strand specific RNAseq reads and UCSC genome
browser custom tracks. *BMC Genomics* **16**, 145
(2015).](http://paperpile.com/b/VjSfVH/ZTg03)

222\. [Wu, E., Nance, T. & Montgomery, S. B. SplicePlot: a utility for
visualizing splicing quantitative trait loci. *Bioinformatics* **30**,
1025--1026 (2014).](http://paperpile.com/b/VjSfVH/F6BgU)

223\. [Rogé, X. & Zhang, X. RNAseqViewer: visualization tool for RNA-Seq
data. *Bioinformatics* **30**, 891--892
(2014).](http://paperpile.com/b/VjSfVH/9Ywn0)

224\. [Tokheim, C., Park, J. W. & Xing, Y. PrimerSeq: Design and
visualization of RT-PCR primers for alternative splicing using RNA-seq
data. *Genomics Proteomics Bioinformatics* **12**, 105--109
(2014).](http://paperpile.com/b/VjSfVH/0lcAx)

225\. [Chelaru, F., Smith, L., Goldstein, N. & Bravo, H. C. Epiviz:
interactive visual analytics for functional genomics data. *Nat.
Methods* **11**, 938--940 (2014).](http://paperpile.com/b/VjSfVH/PdrCV)

226\. [Mariette, J. *et al.* RNAbrowse: RNA-Seq de novo assembly results
browser. *PLoS One* **9**, e96821
(2014).](http://paperpile.com/b/VjSfVH/VoN5Q)

227\. [Severin, J. *et al.* Interactive visualization and analysis of
large-scale sequencing datasets using ZENBU. *Nat. Biotechnol.* **32**,
217--219 (2014).](http://paperpile.com/b/VjSfVH/KL5oc)

228\. [Liu, Q. *et al.* Detection, annotation and visualization of
alternative splicing from RNA-Seq data with SplicingViewer. *Genomics*
**99**, 178--182 (2012).](http://paperpile.com/b/VjSfVH/d5CT4)
