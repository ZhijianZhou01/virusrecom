# VirusRecom: Detecting viral recombination

## 1.Download and install
VirusRecom and all the updated versions is freely available for non-commercial user at https://github.com/ZhijianZhou01/VirusRecom/releases. After obtaining the program,
users could directly run the program in Windows, MacOS or Linux (Ubuntu 16.04 or more) systerms without installation.

In general, the executable file of VirusRecom is located at the Main folder. Then, just double click the VirusRecom.exe (windows system) or VirusRecom (Linux or MacOS system) to start. If you could not get permission to run VirusRecom on linux system, you could change permissions by ```chmod -R 777 VirusRecom```. 


## 2.Getting help
VirusRecom is a command line interface program, users can get help documentation of the software by entering  ```VirusRecom -h ``` or  ```VirusRecom --help ```. 

 ```
usage: 
VirusRecom [-h] [-a ALIGNMENT] [-q QUERY] [-l LINEAGE] [-g GAP] [-m METHOD] 
[-w WINDOW] [-s STEP] [-mr MAX_REGION] [-p PERCENTAGE] [-b BREAKPOINT] 
[-bw BREAKWIN] [-t THREAD] [-y Y_START]
 ```

## 3.
In linux systerms, MAFFT may not work properly，please modify the “prefix path” in mafft program (external_program/mafft/linux/bin/mafft),it might have been so before in file of mafft:

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
```
