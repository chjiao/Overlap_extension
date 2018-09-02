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

## Test data and examples   
The test data sets are in folder test_data/.   
virus.fa   
This data set contains simulated viral reads from HIV-1, HCV genotype 1, and HGV. 
HIV.sam   
This SAM file contains a small subset of aligned HIV-1 reads. With these aligned (287) reads, more HIV-1 reads can be recruited from the viral metagenomics data set.   

Example:    
cd Overlap_extension/   
./build -f test_data/virus.fa -o virus    
./overlap -S test_data/HIV.sam -x virus -f test_data/virus.fa -c 180 -o virus_recruit.fa   
If everything is good, the recruited reads number should be 8008.   
