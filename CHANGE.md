### 1.4.0

Update with a graph theory of breakend, BND=(BP1,BP2), where BP=(chr,pos,dir).

Machine learning modules for ranger and xgboost were replaced with cpp.

Quality control module (fastq) was included in etching_filter

New Pan-Genome K-mer set, PGK2 (based on 894 human genomes), was released (http://big.hanyang.ac.kr/ETCHING/).



### 1.3.7a

Debug to exit if etching_filter stopped abnormally, and to set -f properly in etching (line 793-4)

### 1.3.7

Debug to stop if etching_caller predicted no SV, or etching_sorter removed all SVs.

### 1.3.6c

virtual environment is implemented in etching pipeline to solve dependencies.

### 1.3.6b

etching debug (line 882)
Indentation error fixed (Sorter/scorer_XGBoost)
README updated

### 1.3.6a

File names of final result modified

### v1.3.6

--target-filter and --miscall-kmer-cutoff options were added.


### v1.3.5

Bug fixed (etching line 1283)


### v1.3.4

Debug ```etching``` and ```etching_filter```
They did not run properly when ```-o``` option was not used.


### v1.3.3

```--strand-aware``` and ```--fusion-window``` options added
for fusion-gene detection

### v1.3.2

Updated to stop when an error occurs
README updated
Debug: PGK path misfinding, match_pair

### v1.3.1

Fix option description
etching_fg_identifier update with range


### v1.3.0

Set -M as default
Fix bug for the genome (fa) files with ID lines of multiple words.
Removing artifact step added
Machine learning models updated


### v1.2.3

Calibrating SV length
Accelerating BAM-mode
Debug cut_by_score
Add cut_by_length

