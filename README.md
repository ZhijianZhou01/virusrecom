# VirusRecom: Detecting recombination of viral lineages using information theory

![](https://img.shields.io/badge/System-Windows/Linux/MacOS-green.svg)
![](https://img.shields.io/pypi/pyversions/virusrecom)
[![](https://img.shields.io/badge/Doi-10.1093/bib/bbac513-yellow.svg)](https://doi.org/10.1093/bib/bbac513) 

![](https://img.shields.io/pypi/wheel/virusrecom)
![](https://img.shields.io/pypi/dm/virusrecom)
![](https://img.shields.io/github/stars/ZhijianZhou01/virusrecom)



## 1. Download and install

VirusRecom is developed based on ```Python 3```, and you can get and install VirusRecom in a variety of ways.

### 1.1. pip method

virusrecom has been distributed to the standard library of PyPI, and can be easily installed by the tool ```pip```.

```
pip install virusrecom
virusrecom -h
```

### 1.2. Or local installation

In addition to the  ```pip``` method, you can also install virusrecom manually using the file ```setup.py```. 

Firstly, download this repository, then, run:
```
python setup.py install
virusrecom -h
```

### 1.3. Or run the source code directly

virusrecom can also be run using the source code without installation. First, download this repository, then, install the required python environment of virusrecom:

```
pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple
```

finally, run virusrecom by the file ```main.py```. Please view the help documentation by ```python main.py -h```.

### 1.4. Or use the binary files

For the two earlier release packages (versions v1.0 and v1.1), you can also directly run the binary files of virusrecom without installation. The  binary files are provided at https://github.com/ZhijianZhou01/virusrecom/releases. 

In general, the executable file of virusrecom is located at the  ```main``` folder. Then, running the ```virusrecom.exe``` (windows system) or ```virusrecom``` (Linux or MacOS system) to start. If you could not get permission to run virusrecom on Linux system or MacOS system, you could change permissions by ```chmod -R 775 Directory``` or ```chmod -R 777 Directory```. 


## 2. Getting help
virusrecom is a command-line-interface program, users can get help documentation of the software by entering  ```virusrecom -h ``` or  ```virusrecom --help ```. 

<b>For detailed documentation, please refer to</b> [Manual of VirusRecom v1.1](https://github.com/ZhijianZhou01/virusrecom/blob/main/Manual%20of%20VirusRecom%20v1.1_2023.12.23.pdf)

<b>Tip: since version 1.1, virusrecom optimizes the parameters of input-file, which is slightly different from virusrecom v1.0.</b>

<b>The simple help documentation of virusrecom v1.1 is as follows.</b>

| Parameter | Description |
| --- | --- |
|-h, --help | Show this help message and exit.|
|-a ALIGNMENT | Aligned sequence file (*.fasta). Note, each sequence name requires containing lineage mark.|
|-ua UNALIGNMENT | Unaligned (non-alignment) sequence file (*.fasta). Note, each sequence name requires containing lineage mark.|
|-at ALIGN_TOOL | Program used for multiple sequence alignments (MSA).|
|-iwic INPUT_WIC | Using the already obtained WIC values of reference lineages directly by a *.csv input-file.|
|-q QUERY | Name of query lineage (usually potential recombinant), such as ‘-q xxxx’. Besides, ‘-q auto’ can scan all lineages as potential recombinant in turn.|
|-l LINEAGES | Path of a text-file containing multiple lineage marks.|
|-g GAP | Reserve sites containing gap in subsequent analyses? ‘-g y’means to reserve, and ‘-g n’ means to delete.|
|-m METHOD | Method for scanning. ‘-m p’ means use polymorphic sites only, ‘-m a’ means use all the sites.|
|-w WINDOW | Number of nucleotides sites per sliding window. Note: if the ‘-m p’ has been used, -w refers to the number of polymorphic sites per windows.|
|-s STEP | Step size of the sliding window. Note: if the ‘-m p’ has been used, -s refers to the number of polymorphic sites per jump.|
|-mr MAX_REGION | The maximum allowed recombination region. Note: if the ‘-m p’ method has been used, it refers the maximum number of polymorphic sites contained in a recombinant region.|
|-cp PERCENTAGE | The cutoff threshold of proportion (cp, default: 0.9) used for searching recombination regions when mWIC/EIC >= cp, the maximum value of cp is 1.|
|-cu CUMULATIVE | Simply using the max cumulative WIC of all sites to identify the major parent. Off by default. If required, specify ‘-cu y.|
|-b BREAKPOINT | Possible breakpoint scan of recombination. ‘-b y’ means yes, ‘-b n’ means no. Note: this option only takes effect when ‘-m p’ has been specified.|
| -bw BREAKWIN | The window size (default: 200) used for breakpoint scan. The step size is fixed at 1. Note: this option only takes effect when ‘-m p -b y’ has been specified.|
|-t THREAD | Number of threads (default: 1) used for MAS.|
|-y Y_START | Starting value (default: 0) of the Y-axis in plot diagram.|
|-le LEGEND | The location of the legend, the default is adaptive. '-le r' indicates placed on the right.|
|-owic ONLY_WIC | Only calculate site WIC value. Off by default. If required, please specify ‘-owic y’.|
|-e ENGRAVE | Engraves file name to sequence names in batches. By specifying a directory containing one or multiple sequence files (*.fasta).|
|-en EXPORT_NAME | Export all sequence name of a *.fasta file.|
|-o | Output directory to store all results.|


For more information about the algorithm of virusrecom, please refer to [the publication of virusrecom](https://academic.oup.com/bib/article-abstract/24/1/bbac513/6886420).

## 3. Example of usage
The sequences data for test in the documentation was stored at https://github.com/ZhijianZhou01/virusrecom/tree/main/example. 

<b>Note, the ```recombination_test_data.zip``` in directory ```example``` is against virusrecom v1.0, not virusrecom v1.1</b>.

In this demonstration, the test data is from the the ```recombination_test_data_v1.1.zip``` provided in the directory ```example```. 

### 3.1. Aligned input-sequences
If the input sequence-data has been aligned, and it should be loaded via the ```-a``` parameter. Multiple sequence alignments (MSA) can be pre-completed by many programs, this is not introduced. Now, let's focus on the directory ```aligned_input_sequences``` in the file ```recombination_test_data_v1.1.zip```. 

(1) An aligned sequence-file named ```alignment_lineages_data.fasta```, which including multiple sequences from the query lineage and other reference lineages. 
    
(2) A text-file named ```reference_lineages_name.txt```, which including the names (marks) of these reference lineages. 
    
```
     reference_lineage_1
     reference_lineage_2
     reference_lineage_3
     reference_lineage_4
     reference_lineage_5
     reference_lineage_6
     reference_lineage_7
     reference_lineage_8
     reference_lineage_9
```
     
Note, these marks of reference lineages should also appear in sequence names of the file ```alignment_lineages_data.fasta```. <b>The mark of each reference lineage should be unique</b>, otherwise, there will be duplicate matches in subsequent analysis.

Before running the command of VirusRecom, let's think about the search strategy for recombination events. Firstly, we use only polymorphic sites considering that sequences from these lineages are highly similar, which means that the parameter ```-m p``` needs to be specified. Secondly, we do not consider gap-containing sites in this test and use the parameter ```-g n```. Instead, if you consider these gap sites, you need to use the parameter ```-g y```. Next, in the first run, let's try first with a window size of 100 and a step size of 20. Of note the value of “size” at this time represents the number of polymorphic sites because the ```-m p``` parameter has been specified. For the two parameters ```-cp``` and ```-mr```, we use the default value of 0.9 and 1000 in this test. Finally, we specify a folder to save the results by parameter ```-o```. 

Then, switch the current directory to ```aligned_input_sequences```, and run the following command (an example) to detect recombination events in query lineage:

```
virusrecom -a alignment_lineages_data.fasta -q query_recombinant -l reference_lineages_name.txt -g n -m p -w 100 -s 20 -o outdir
```

Note: (1) if the current directory is not switched to ```aligned_input_sequences```, the file and directory path in command need the absolute paths instead of relative paths.
(2) the string “query_recombinant” in command is the corresponding mark of query lineage in the file ```alignment_lineages_data.fasta```.


<b>After the run is complete</b>, in the directory ```outdir```, there are three subdirectories and two aggregated reports:

![outdir.png](https://github.com/ZhijianZhou01/virusrecom/blob/main/figture/outdir.png)

(1) In the directory ```run_record```, if ```-g n``` is specified, and the file ```Record_of_deleted_gap_sites_*.txt``` containing all the gap sites will be created. Besides, If ```-m p``` is specified, and the file ```Record_of_same_sites_in_aligned_sequence*.txt``` containing all the same sites will be created.

(2) In the directory ```WICs_of_sites```, the file ```*_site_WIC_from_lineages.pdf```, ```*_site_WIC_from_lineages.xlsx``` and the file ```*_site_WIC.csv``` are used to record the WIC value of each site. 

(3) In the directory ```WICs_of_slide_window```, the file ```*_mWIC_from_lineages.xlsx``` and the file ```*_mWIC_from_lineages.pdf``` are used to record the mean WIC of each sliding window. 

![recombination_step4.png](https://github.com/ZhijianZhou01/virusrecom/blob/main/figture/recombination_step4.png)

The user can fine-tune the window size and step size according to the density of points in the generated graph. In general, very dense points means that the noise is too high and the window size can be increased appropriately in next scan. 

In addition to the three sub-directories above, VirusRecom provides two summary files. The file ```Possible_recombination_event_conciseness.txt``` only retains results of recombination events with p-values less than 0.05.
```
Possible major parent: reference_lineage_1(global mWIC: 1.8976186779157704)

Other possible parents and significant recombination regions (p<0.05):
reference_lineage_2	7237 to 11539(mWIC: 1.9553354371515168), p_value: 7.831109305531908e-06	

Significance test of recombinant regions using Mann-Whitney-U test with two-tailed probabilities, p-value less than 0.05 indicates a significant difference.
```

In this output report, the major parent of query lineage was ```reference_lineage_1``` and the minor parent was ```reference_lineage_2```, and the recombination region was site 7237 to 11539 and the p-value was 7.83e-06. The identified recombination event was relatively close to the actual (from site 7333 to 11473 in the genome), and the error of the recombination boundary is also acceptable.

In fact, ```Possible_recombination_event_conciseness.txt``` is interpretations of the recombination information contained in ```*_mWIC_from_lineages.pdf```. Although VirusRecom shows a good balance between precision and recall in simulated data, false positive or false negatives sometimes occur. Therefore, for the identification results from VirusRecom, users can make own judgment. 

Besides, the output file ```Possible_recombination_event_detailed.txt``` shows those results with p-values greater than 0.05. <b>Tip: recombination events with p-values over 0.001 are less reliable</b>. 

If ```-b y``` is specified, then VirusRecom will perform the search of recombination breakpoint and plot. For example:
```
virusrecom -a alignment_lineages_data.fasta -q query_recombinant -l reference_lineages_name.txt -g n -m p -w 100 -s 20 -b y -bw 200 -o outdir
```
<b>Tip:</b> (1) ```-b y``` only takes effect when ```-m p``` has been specified. 
(2) the step size of breakpoint search is fixed to 1. 

The negative logarithm of p-value in each site is in the file ```*_-lg(p-value)_for_potential_breakpoint.pdf``` and the file ```*_-lg(p-value)_for_potential_breakpoint.xlsx```. 

![breakpoint.jpg](https://github.com/ZhijianZhou01/virusrecom/blob/main/figture/breakpoint.jpg)

The highest peak (the highest −lgP value) indicated the possible recombination breakpoint.

### 3.2. Unaligned input-sequences
VirusRecom can also handle unaligned input-sequences. In this case, multiple sequence alignment is performed by calling external program. In <b>virusrecom v1.1</b>, mafft, muscle, and clustal-omega is supported. It is worth mentioning that VirusRecom call them from the system path, so they need to be installed on the machine beforehand.

For the example data in directory ```unaligned_input_sequences```, run the following command:
```
virusrecom -ua unalignment_lineages_data.fas -at mafft -q query_recombinant -l reference_lineages_name.txt -g n -m p -w 100 -s 20 -o outdir
```
<b>Note:</b> (1) ```-at mafft``` means to call mafft in the system path, and the alignment strategy is auto. Besides, using ```-at muscle``` to call muscle and using ```-at clustalo``` to call clustal-omega.
(2) the string ```query_recombinant``` in command is the corresponding mark of query lineage in the file ```unalignment_lineages_data.fas```.

The interpretation of the output result is consistent with section 3.1. 


### 3.3. Non-lineage data
In VirusRecom, the reference lineage is allowed to contain only one single sequence. Under this condition, mWIC value of the fragment is essentially a multiple of shared identity. If -g n is used in the calculation, the mWIC is twice as large as shared identity. If ```-g y``` is used in the calculation, the mWIC is $\log_2{5}$ as large as shared identity. 

Of noted, for recombination analysis without lineage data, the additional feature is only recommended for non-highly similar sequences and the user can use it to draw an identity point map.

The test data is in directory ```non_lineage_data``` of the file ```recombination_test_data_v1.1.zip```. 

The Delta-CoV HNU1-1 is a known recombinant from SpCoV HKU17-USA and ThCoV HKU12, and the break points were identified at genome positions nt 21017 and 25056, which is jointly identified and confirmed by RDP3 and Simplot by [Wang et al., 2022](https://onlinelibrary.wiley.com/doi/10.1111/tbed.14029). 

Considering that they are not highly similar sequences, we use all sites (```-m a```) in the alignment. Then, we use a larger window value, and run following command:
```
virusrecom -a alns.fasta -q HNU1-1 -l alns_seq_taxon.txt -g n -m a -w 800 -s 100 -cp 0.7 -mr 6000 -le r -o output

```

The mWIC from reference lineages is as follows:

![hnu1-1.jpg](https://github.com/ZhijianZhou01/virusrecom/blob/main/figture/hnu1-1.jpg)

<b>Note,</b> because each “lineage” contains only one sequence and ```-g n``` is used in the example, the mWIC in the picture is actually twice the size of “sequence identity”. 

The possible recombination event identified by VirusRecom is as follows:
```
Possible major parent: HKU17-USA(global mWIC: 1.5914816042426252)  

Other possible parents and significant recombination regions (p<0.05): 
HKU12	20720 to 25297(mWIC: 1.8039433490697028), p_value: 2.783880536189705e-204
```

The possible major parent of HNU1-1 is HKU17-USA and minor parent is HKU12, and the recombination region is about 20720-25297 nt in the alignment.


## 4. Common questions
### 4.1. Default values of parameter 
For the value of a parameter, if not specified, the software uses the default value. 
However, the default value is not suitable for all data. In addition to window size (```-w```) and step size (```-s```) of sliding window, values of ```-cp``` and ```-mr``` also require users to adjust based on the data. 

When VirusRecom runs, the value of each parameter is printed printed on the screen and you can check them. 

### 4.2. How to mark lineage in sequence name?
Typically, this is part of the data preparation. In virusrecom v1.1, users can easily get it done via ```-e``` parameter. The ```-e``` parameter can engrave file-name to sequence names in batches. The example is as follows:
```
virusrecom -e input_directory -o outdir
```
<b>Tip:</b> The directory ```input_directory``` can contain multiple fasta files, and each fasta file can contain multiple sequences. After the running, finally, each sequence name will contain its file-name. 

Therefore, if the file-name of fasta file is a lineage name, the lineage name can be written into the sequence name in batches.

### 4.3. How to change the color scheme in an image?
If you own programming skills, you can directly modify the order of the colors in the ```plt_corlor_list.py``` file. If not, you can use output matrix provided by VirusRecom, and they are usually suffixed with ```.xlsx```. 


## 5. Attention
(1) If installed software was from <b>the binary files of virusrecom v1.0</b> (https://github.com/ZhijianZhou01/virusrecom/releases/tag/v1.0) , and need to call the plug-in MAFFT for multiple sequence alignment in the linux systerms, MAFFT may not work properly，please modify the “prefix path” in mafft program (external_program/mafft/linux/bin/mafft),it might have been so before in file of mafft:

```
if [ "$MAFFT_BINARIES" ]; then
	prefix="$MAFFT_BINARIES"
else        
	prefix=../external_program/mafft/linux/libexec/mafft
fi
```

Then you need to modify it to something like this (for example, directory of "virusrecom_V1.0_linux" was directly placed in the ```/home```):

```
if [ "$MAFFT_BINARIES" ]; then
	prefix="$MAFFT_BINARIES"
else        
	prefix= /home/VirusRecom_V1.0_linux/external_program/mafft/linux/libexec/mafft
fi
```

(2) If installed software was from <b>the binary files of virusrecom v1.1</b> (https://github.com/ZhijianZhou01/virusrecom/releases/tag/v1.1), or from some other form (including ```pip install virusrecom``` and run the source code of virusrecom directly), when the unaligned input-sequences was used in analysis, the external program (such as MATTF, MUSCLE) used for multiple sequence alignment need to be installed and added to environment variables of system or user beforehand, because VirusRecom call them from the environment variables directly.


## 6. Citation
Zhou ZJ, Yang CH, Ye SB, Yu XW, Qiu Y, Ge XY. VirusRecom: an information-theory-based method for recombination detection of viral lineages and its application on SARS-CoV-2. <i>Brief Bioinform</i>. 2023 Jan 19;24(1):bbac513. [doi: 10.1093/bib/bbac513](https://academic.oup.com/bib/article-abstract/24/1/bbac513/6886420). PMID: 36567622.
