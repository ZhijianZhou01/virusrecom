# VirusRecom: Detecting recombination of viral lineages using information theory

## 1. Download and install
Users can run the source code directly (such as ```python main.py -h```), or use the releases we provided at https://github.com/ZhijianZhou01/VirusRecom/releases. VirusRecom and all the updated versions is freely available for non-commercial. After obtaining the program, users could directly run the program in Windows, MacOS or Linux systerms without installation.

Take the releases as an example, in general, the executable file of VirusRecom is located at the  ```main``` folder. Then, running the VirusRecom.exe (windows system) or virusrecom (Linux or MacOS system) to start. If you could not get permission to run VirusRecom on Linux system or MacOS system, you could change permissions by ```chmod -R 777 Directory```. 


## 2. Getting help
VirusRecom is a command line interface program, users can get help documentation of the software by entering  ```VirusRecom -h ``` or  ```VirusRecom --help ```. 

 ```
usage: virusrecom [-h] [-a ALIGNMENT] [-q QUERY] [-l LINEAGE] [-g GAP]
                  [-m METHOD] [-w WINDOW] [-s STEP] [-mr MAX_REGION]
                  [-cp PERCENTAGE] [-b BREAKPOINT] [-bw BREAKWIN] [-o OUTDIR]
                  [-t THREAD] [-y Y_START]

optional arguments:
  -h, --help      show this help message and exit
  -a ALIGNMENT    FilePath of an aligned sequence set(*.fasta format)
                  containing all sequences used for analysis, then the
                  alignment will be skipped. Default is null. If using, name
                  of each sequence in aligned sequence set requires containing
                  the mark(a unique string) of the lineage.
  -q QUERY        FilePath of query lineage (potential recombinant, *.fasta
                  format). Note, if the '-a alignment' has been used, please
                  enter the mark (a unique string) of queried recombinant
                  here, such as '-q XE_'(not a FilePath), besides, using '-q
                  auto' and all sequences will be scanned as potential
                  recombinants in turn.
  -l LINEAGE      DirPath of reference lineages. One sequence file (*.fasta
                  format) per lineage, and each lineage could contain multiple
                  sequences. Note, if the '-a alignment' has been used, please
                  enter a text file containing the marks (a unique string) of
                  lineages here, not a DirPath.
  -g GAP          Gaps (-) in the alignment were used in analysis? '-g y':
                  reserve gaps, '-g n': delete gaps.
  -m METHOD       Scanning method of recombination analysis. '-m p': using
                  polymorphic sites only, '-m a': using all the monomorphic
                  sites and polymorphic sites.
  -w WINDOW       Number of nt sites per sliding window. Note: if the '-m p'
                  method has been used, -w refers to the number of polymorphic
                  sites per windows.
  -s STEP         Step size for scanning these sites. Note: if the '-m p'
                  method has been used, -w refers to the number of polymorphic
                  sites per jump.
  -mr MAX_REGION  The maximum allowed recombination region. Note: if the '-m
                  p' method has been used, it refers the maximum number of
                  polymorphic sites contained in a recombinant region.
  -cp PERCENTAGE  The cutoff threshold of proportion (cp, default was 0.9) for
                  searching recombination regions when mWIC/EIC >= cp, the
                  maximum value of cp is 1. For detection in genus level,
                  about 0.5 is recommended.
  -b BREAKPOINT   Whether to run the breakpoint scan of recombination. ‘-b y’:
                  yes, ‘-b n’: no. Note: this option only takes effect when
                  '-m p' has been specified!
  -bw BREAKWIN    The window size (polymorphic sites, default is 200) used for
                  breakpoint scan. The step size is fixed at 1. Note: this
                  option only takes effect when '-m p -b y' has been
                  specified!
  -o OUTDIR       The outdir of results.
  -t THREAD       Number of threads used for the multiple sequence alignments
                  (MSA), default is 1.
  -y Y_START      Specify the starting value of the Y axis in the picture, the
                  default is 0.

-----------------------------------------------------
☆ Example of use ☆
  (1) If the input-sequence data was not aligned:
      virusrecom -q XE.fasta -l Lineage_Dir -g n -m p -w 100 -s 20 -t 2

  (2) If the input-sequence has been aligned:
      virusrecom -a alignment.fasta -q XE_ -l lineage_name_list.txt -g n -m p -w 100 -s 20
 ```

## 3. Attention
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
