# VirusRecom: Detecting recombination of viral lineages using information theory

## 1. Download and install
Users can run the source code directly (such as ```python main.py -h```), or use the releases we provided at https://github.com/ZhijianZhou01/VirusRecom/releases. VirusRecom and all the updated versions is freely available for non-commercial. After obtaining the program, users could directly run the program in Windows, MacOS or Linux systerms without installation.

Take the releases as an example, in general, the executable file of VirusRecom is located at the  ```main``` folder. Then, running the VirusRecom.exe (windows system) or virusrecom (Linux or MacOS system) to start. If you could not get permission to run VirusRecom on Linux system or MacOS system, you could change permissions by ```chmod -R 777 Directory```. 


## 2. Getting help
virusrecom is a command line interface program, users can get help documentation of the software by entering  ```virusrecom -h ``` or  ```virusrecom --help ```. 

 ```
-h, --help	Show this help message and exit.
-a	FilePath of an aligned sequence set (*.fasta format) containing all sequences used for analysis, then the sequence alignment will be skipped. Default value is null. If “-a” parameter was used, the name of each sequence in aligned sequence set requires containing the mark (a unique string) of the lineage.
-q	FilePath of query lineage (usually potential recombinant, *.fasta format). Note, if the ‘-a’ parameter has been used, please enter the mark (a unique string) of query lineage here, such as ‘-q xxxx’, not a FilePath. Using ‘-q auto’ and all lineages will be scanned as potential recombinants in turn.
-l	DirPath of reference lineages. One sequence file (*.fasta format) per lineage, and each lineage could contain multiple sequences. Note, if the ‘-a’ parameter has been used, please enter a file path of a text file containing the mark (a unique string) of lineage here, not a DirPath.
-g	Gaps (-) in the alignment were used in analysis? ‘-g y’ means to reserve gaps, and ‘-g n’ means to delete gaps.
-m	Scanning method of recombination analysis. ‘-m p’ means use polymorphic sites only, ‘-m a’ means all the monomorphic sites and polymorphic sites.
-w	Number of nucleotides sites per sliding window. Note: if the ‘-m p’ has been used, -w refers to the number of polymorphic sites per windows.
-s	Step size of the sliding window. Note: if the ‘-m p’ has been used, -s refers to the number of polymorphic sites per jump.
-mr	The maximum allowed recombination region. Note: if the ‘-m p’ method has been used, it refers the maximum number of polymorphic sites contained in a recombinant region.
-cp	The cutoff threshold of proportion (cp, default value was 0.9) used for searching recombination regions when mWIC/EIC >= cp, the maximum value of cp is 1.
-cm	Whether to simply use the max cumulative WIC of all sites to identified the major parent. The default value is ‘n’ and means ‘no’. If required, please specify ‘-cm y’.
-b	Whether to run the breakpoint scan of recombination. ‘-b y’ means yes, ‘-b n’ means no. Note: this option only takes effect when ‘-m p’ has been specified.
-bw	The window size (polymorphic sites, default value is 200) used for breakpoint scan. The step size is fixed at 1. Note: this option only takes effect when ‘-m p -b y’ has been specified.
-t	Number of threads used for the multiple sequence alignments (MSA), default is 1.
-y	Specify the starting value of the Y-axis scale in the picture, the default is 0.

```

## 3. Example of usage
The sequences data for test in the manual was stored at https://github.com/ZhijianZhou01/virusrecom/tree/main/example. Take the ```recombination_test_data.zip``` provided in the directory ``example`` as a demonstration.

### 3.1. Unaligned input-sequences
VirusRecom owns a pipeline built in to handle unaligned sequences. In this case, multiple sequence alignment is performed by MAFFT (Katoh & Standley, 2013) with alignment strategy of “auto” parameter. In the directory of unaligned_input_sequences in the compressed file of recombination_test_data.zip, and the query lineage (simulated recombinant) including multiple sequences was in the file query_recombinant.fasta, and these reference lineages were placed in a separate directory named lineages_dir. For the reference lineages, sequences for each lineage need to be placed in a separate file, such as reference_lineage_1.fasta, reference_lineage_2.fasta and reference_lineage_3.fasta in the lineages_dir directory. Of note, the name of each file is an important label to distinguish different lineages. In fact, the query lineage in the test data was known recombinant because it was from the synthetic data, the major parent was reference_lineage_1, the minor parent was reference_lineage_2 and the recombination region was from site 7333 to 11473 in the genome. However, here we treat the query lineage as a potential recombinant.















## 4. Attention
If you need to call MAFFT for multiple sequence alignmentIn in linux systerms, MAFFT may not work properly，please modify the “prefix path” in mafft program (external_program/mafft/linux/bin/mafft),it might have been so before in file of mafft:

```
if [ "$MAFFT_BINARIES" ]; then
	prefix="$MAFFT_BINARIES"
else        
	prefix=../external_program/mafft/linux/libexec/mafft
fi
```

Then you need to modify it to something like this (for example, directory of "VirusRecom_V1.0_linux" was directly placed in the ```/home```):

```
if [ "$MAFFT_BINARIES" ]; then
	prefix="$MAFFT_BINARIES"
else        
	prefix= /home/VirusRecom_V1.0_linux/external_program/mafft/linux/libexec/mafft
fi
```
