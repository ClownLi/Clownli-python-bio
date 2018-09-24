#GATK是Genome Analysis Toolkit的缩写，是用来处理高通量测序数据的一套软件。
#最初，GATK被设计用来分析人类基因组和外显子，主要用来寻找SNP和indel。
#后开，GATK的功能越来越丰富，增加了short variant calling、计算copy number（CNV）和结构变异（SV）等新功能。
#同时，GATK也越来越广泛地应用于其他物种的数据分析中。现在，GATK已经成为了基因组和RNA-seq分析过程中，寻找变异的行业标准。

#用GATK寻找SNP和Indel，有一个标准的分析流程叫做GATK Best Practise主要包括以下几个步骤：
##1 原始数据质控 （原始测序数据 ---> QC/过滤低质量read数据）
##2 数据预处理   （read比对 ---> 排序(sort) ---> 去重复 ---> 局部重比对 ---> 碱基质量重校正(BQSR)）
##3 变异检测    （变异检测 ---> Merge(optional) ---> (joint Genotype)
#              （变异检测 ---> 变异质控和过滤（VQSR）                        ---> 最终变异结果）

####http://www.huangshujia.me//2017/09/19/2017-08-25-Begining-WGS-Data-Analysis-The-pipeline.html


#1 对原始下机fastq文件进行过滤和比对
#a 对参考基因建立索引
bwa index ref.fa
#b BWA Alignment
bwa mem -t 24 -R "@RG\tID:<ID>\tLB:<LIBRARY_78>\tSM:<78>\tPL:ILLUMINA" ref.fa read1.fq read2.fq > sample.sam

#GATK对基因组要求一个字典文件
java -jar /home/stu_madongna/software/picard.jar CreateSequenceDictionary R=ref.fasta O=ref.dict

#samtools faidx 能够对fasta 序列建立一个后缀为.fai 的文件，根据这个.fai 文件和原始的fastsa文件，能够快速的提取任意区域的序列
samtools faidx /home/stu_madongna/SOAPDATA/ref.fa

#2 sam文件转bam文件
samtools view -bS sample.sam -o sample.bam

#3 对bam文件排序
java -jar /home/lina/miniconda3/share/picard-2.18.14-0/picard.jar SortSam I=sample.bam O=sample.sort.bam SORT_ORDER=coordinate

#4 对bam文件进行加头(head)处理
java -jar /home/lina/miniconda3/share/picard-2.18.14-0/picard.jar AddOrReplaceReadGroups I=sample.sort.bam O=sample.sort.addhead.bam ID=sampleID LB=sampleID PL=illumina PU=samplePU SM=sample

#5 标记重复序列
java  -Xmx15g -jar /home/lina/miniconda3/share/picard-2.18.14-0/picard.jar MarkDuplicates I=sample.sort.addhead.bam O=sample.rmdup.bam REMOVE_DUPLICATES=false  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 METRICS_FILE=sample.sort.addhead.bam.metrics

#6 要对上一步得到的结果生成索引文件
samtools index sample.rmdup.bam

#7 Local realignment around indels
#a 将比对到indel附近的reads进行局部重新比对，将比对的错误率降到最低
java -Xmx50g -jar /home/stu_madongna/software/GenomeAnalysisTK.jar -R ref.fa -T RealignerTargetCreator -I sample.rmdup.bam -o sample.realign.intervals
#b 通过运行IndelRealigner在这些区域内进行重新比对
java -Xmx50g -jar /home/stu_madongna/software/GenomeAnalysisTK.jar -R ref.fasta -T IndelRealigner -targetIntervals sample.realign.intervals -I sample.rmdup.bam -o sample.realign.bam

#9 GATK的HaplotypeCaller进行Variant calling   ### 如果做群体SNP 在-I后面加多个bam文件
java -Xmx50g -jar /home/stu_madongna/software/GenomeAnalysisTK.jar -R ref.fasta -T HaplotypeCaller -I sample.realign.bam -o sample.gatk.raw.vcf -stand_call_conf 30

#1（另外一种变异检测的方法）使用Samtools检测变异
samtools index sample.realign.bam

#2 call SNP           ### 如果做群体SNP 在-I后面加多个bam文件
samtools mpileup  -s -t DP -t AD -t ADF -t ADR -t INFO/AD -t SP -uvmFOB -C 50 -f  sample.realign.bam | bcftools view --types snps,indels - > sample.samtools.raw.vcf

#### 过滤
perl /home/stu_madongna/perl/VCF_cleaner.pl
