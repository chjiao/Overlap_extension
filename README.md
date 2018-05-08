# Overlap_extension

May 8, 2018 updates:   
1. Allow recruiting paired-end reads at the end   
2. Enable the user to define the number of partitions for indexing   

The codes requires C++11 for compilation. 

## Installation
For compile:     
cd Overlap_extension   
make    

## Usage
After compilation, there will be two binary files: build and overlap   
To run the program:   
(1) build reads index   
./build -f reads.fa -o prefix   

The default number of partitions is 5   

(2) recruite reads    
./overlap -S align.sam -x prefix -f reads.fa -c overlap_cutoff -o recruited_reads.fa    
align.sam is the alignment results of reads.fa on available reference   
