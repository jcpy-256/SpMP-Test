# SpMp-Test 
This project is forked from [SpMP](https://github.com/jcpy-256/SpMP) and we rectify the testing code for SpTRSV, ILU0, IC0. 

## Preliminaries
* GCC >= 8.5 
* OpenMP 

## Getting Started
**1. Build**
We provide convenient scripts (`./script` dir) to handle the build process:
```shell 
# Compile the project
./ make.sh
```
**2. Data Preparation**
This code expects matrix files to be located in the /data directory.

* Place your .mtx files in `/data`.

* Update `matrixList.txt` to include the filenames of the matrices you wish to test.

**4. Running**
We have provided specialized scripts for different sparse solvers. You can run them directly from the `script` directory.
```shell
./run_trsv.sh # SpTRSV 
./run_ic.sh   # SpIC0
./run_ilu.sh  # SpILU0
```
**Note**: You can customize the **thread count** and **iteration rounds** within these scripts to match your hardware.

**5. Output & Results:**
All generated files are stored in different directories:

`csv/`: Contains the measured performance metrics and benchmark results.

`log/`: Contains detailed program execution logs.

`out/`: Contains standard output redirected from the running scripts.





