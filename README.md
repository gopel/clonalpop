Docker image of ClonalPop, an automatic pipeline for analysis of bacterial clonal populations created by Yanis Bendjelal  
GitHub :  https://github.com/gopel

# Versions 

ClonalPop: `1.0`  
SPAdes :  `3.13.0`    
Unicycler : `0.4.7`  
Prokka : `1.13`  
Breseq : `0.33.2`    

# Installation

## With Docker (the easiest way) 

### Download and install Docker on your machine :  
https://www.docker.com/products/docker-desktop  

### Download this image (copy/paste in your Terminal) :  
```
docker pull gopel/clonalpop  
```

## Without docker
brew

# Usage 
```
docker run gopel/clonalpop path/to/input/file path/to/output/file
```   

At the end of the pipeline, all the intermediary files are automatically supressed. To keep one, add :  
SPAdes :  `spades`    
Unicycler : `unicycler`    
Prokka : `prokka`    
Breseq : `breseq`      

For example, to keep comparison files :   
`docker run gopel/clonalpop path/to/input/file path/to/output/file breseq`

# Hijack the pipeline

### Assembly pipeline

To assemble Illumina reads :    
```
docker run gopel/clonalpop path/to/input/file path/to/output/file spades
``` 

To assemble Nanopore files with associated Illumina reads :    
```
docker run gopel/clonalpop path/to/input/file path/to/output/file unicycler
```  

### Annotation pipeline

Give the files you want to annotate with the extension '_reference', then enter :    
 ```
docker run gopel/clonalpop path/to/input/file path/to/output/file prokka
```

### Characterisation pipeline

To characterise a file' reference content :  
 ```
docker run gopel/clonalpop path/to/input/file path/to/output/file abricate
```

### Comparison pipeline

To keep the comparison files :   
 ```
docker run gopel/clonalpop path/to/input/file path/to/output/file breseq
```
 
 
