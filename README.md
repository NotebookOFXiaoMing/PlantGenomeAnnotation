# PlantGenomeAnnotation


需要注释的基因组放到一个文件夹下，用到的转录组数据放到另外一个文件加下

step01.yaml 里用绝对路径写基因组的存放位置，转录组的存放位置

genomeAnnotationStep01.smk 做的是 fastp hisat2 samtools stringtie transdecoder 流程生成bam文件用于braker2去训练augustus和genemark transdecoder用于生成用于evm的转录本证据

运行这一步需要用到的软件都可以用conda安装

新建一个叫rnaseq的环境

```
conda create -n rnaseq
```

安装软件

```
conda activate rnaseq
conda install fastp
conda install hisat2
conda install samtools=1.15 # samtools 安装最新版使用的时候总是报错，所以安装的时候指定一个相对比较老的版本
conda install stringtie
conda install transdecoder
conda install snakemake
```

运行第一步

```
snakemake -s genomeAnnotationStep01.smk --configfiles=step01.yaml --cores 128 -p
```

genomeAnnotationStep02.smk 是用上一步转录组生成的bam文件去做braker2

这一步需要用到braker2软件，我单独建了一个环境

```
conda create -n braker2
```

然后用mamba去安装braker2

```
conda activate braker2
conda install mamba
mamba install braker2
mamba install snakemake
```

运行braker2需要配置genemark 

参考一下这个链接 https://www.jianshu.com/p/047216b469e7

运行命令

```
snakemake -s genomeAnnotationStep02.smk --configfiles=step02.yaml --cores 128 -p
```

genomeAnnotationStep03.smk 是使用miniprot将蛋白证据比对到基因组，蛋白证据放到同一个文件夹下，miniprot可以直接使用conda安装，我是安装到了第一步rnaseq的环境下

```
conda activate rnaseq
conda install miniprot
```

运行命令

```
snakemake -s genomeAnnotationStep03.smk --configfile=step03.yaml --cores 48 -p
```

genomeAnnotationStep03_04.smk 是用来生成EVM需要用到的权重文件的，这一步好好想想，应该直接放到step04里

运行命令

```
snakemake -s genomeAnnotationStep03_04.smk --configfile=step04.yaml --cores 48 -p
```

genomeAnnotationStep04.smk 就是evm整合各种证据了，evm也可以直接用conda安装，我是新建了EVM环境

```
conda create -n EVM
conda activate EVM
conda install mamba
mamba install evm
mamba install snakemake
```

运行命令

```
snakemake -s genomeAnnotationStep04.smk --configfile=step04.yaml --cores 128 -p
```


以上内容为运行测试，所以没有准备集群提交任务的脚本，采用的是交互运行的方式进行的，集群提交任务的脚本后续补上
