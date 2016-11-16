# fastARG
Fast heuristic construction of a ancestral recombination graphs from sequence data. 

##Installation
Simply run `make` within the fastARG directory. This requires a C compiler and the zlib library.


##Basic usage

To construct an ARG from haplotype data

	$ fastarg build ex1.txt > output.arg

To output haplotypes from this ARG

	$ fastarg leafseq output.arg > output.txt

Other possible commands are displayed by running `./fastARG` with no additional arguments.

##Input format

For _m_ loci (sites) and _n_ haplotypes/genotypes, the input file should consist of _m_ lines, each giving an (integer) position followed by whitespace, then a string of _n_ non-whitespace characters denoting the state of the _n_ genotypes (or haplotypes) at that locus.

The file `ex1.txt` in the fastARG directory gives an example for _m_=3 loci and _n_=4 haplotypes:

```
1	1000
2	0011
3	0101
```

Note that the first integer on the line (the "position") is currently treated only as a label, and otherwise ignored. The position along the genome is given by the position of the line in the file: the first line denotes the state of haplotypes at position 0, the second the stare at position 1, etc. 

The sequence data consists of non-whitespace characters with character states represented as `0` or `1`. For unphased (genotype) data additionally the character `2` indicates a heterozygote - if no `2`s are present, the characters are assumed to come from haploid genomes. Any other character is taken as missing data (internally represented as `3`).


#Output format

Each line of the ouput file starts with a single letter which indicates the interpretation of the (tab-separated) digits listed on the rest of that line. The initial lines give the starting parameters:

###Starting parameters

###### `E` sEed for random number generator
A single line starting with `E` is followed by an single (long) integer giving the random seed initially passed to the `srand48()` function
   
######  `N` Numerical summary of parameters
A single line starting with `N` is followed by two integers giving the number of loci (_m_) and the number of haplotypes (_n_) respectively

###ARG specification

The following set of lines describe the ARG. Note that nodes in the ARG are indexed from 0...(_n_-1) for leaf nodes, with numbers from _n_ upwards used for internal ARG nodes.

###### `C` Coalescence event
Any number of lines starting with `C` may exist, which indicate branches that coalescence. The initial `C` is followed by 5 or more tab-delimited integers. The first 4 integers represent the current node ID, a descendant node ID, and the half-closed genomic interval over which this coalescence applies. In the normal case that 2 lineages coalesce, there will be a pair of coalescence lines, such as in the following example, which indicates that node 5 consists of a coalescence between nodes 4 and 3

```
C	5	4	0	2	0
C	5	3	0	3	0
```

Note that not all the sequence will coalesce at this point. The parts of the haplotype that coalesce are indicated by the 3rd and 4th number. In this example, the loci that coalesce from node 4 into node 5 are indexed by the interval [0,2) i.e. loci 0 and 1, whereas the parts of the haplotype that coalesce from node 3 into node 5 are in the interval [0,3) i.e. loci 0, 1, and 2.

The fifth number on the coalescence line denotes the number of mutations that are inferred to have happened along the branch connecting this node to the specified descendant. The locus of each mutation is then given by the subsequen numbers. In other words, if the 5th number is 0, it will not be followed by any other numbers (no mutations on this edge). If the 5th number is 1, it will be followed by a single number specifying the locus at which a mutation occurs. If the 5th number is 2, it will be followed by two numbers giving such indices, and so on. 

###### `R` Recombination event
Any number of lines starting with R may exist, which indicate recombination events. This is followed by 5 or more numbers interpreted in the same way as for a coalescence event.

###Sequence specification

###### `S` Sequence at root

A single line starting with an `S` is followed by the ID of the deepest node in all trees (the root), then a sequence of _m_ characters giving the inferred haplotype at this root. This serves to completely define the ancestral locus history.

#Example

The `ex1.txt` file should produce something like the following output

```
E	11
N	4	3
R	4	2	0	3	0
C	5	4	0	2	0
C	5	3	0	3	0
C	6	0	0	3	1	0
C	6	4	2	3	0
C	7	1	0	3	0
C	7	5	0	3	1	1
C	8	6	0	3	0
C	8	7	0	3	1	2
S	8	000
```

showing that the random seed is 11 (line 1) and the number of haplotypes and number of loci is 4 and 3 respectively (line 2).

The ARG indicated by the `R` and `C` records in the above file is shown by the diagram to the left below. Recombination nodes have been prefixed with R, coalescence events with C, and dots represent lines that visually cross, but do not intersect. The diagram to the right is the same but with mutations overlain: `*0` represents a mutation at locus 0, `*1` at locus 1, and `*2` at locus 2

```
       C8                  C8    
      /  \                /  *2   
     /   C7              /   C7  
    /   /  \            /   /  *1 
   C6  /   C5          C6  /   C5
  /  \.   / |         /  \.   / |
 /   .\  /  |        /   .\  /  |
/   /  4R   |       *0  /  4R   |
|   |   |   |       |   |   |   |
0   1   2   3       0   1   2   3

  == ARG ==        == ARG with ==
                   == mutations==
```
The output also states that the line joining node 5 to node 4 only applies to loci in the range [0,2), and the line joining node 4 to node 6 applies to loci in the range [2,3). Thus loci 0 and 1 are inferred to be related as shown in the left hand coalescence tree below, whereas the inferred tree for locus 2 is the right hand tree below.

```  
       C8                C8     
      /  \              /  *2    
     /   C7            /   C7   
    /   /  *1         /   /  \  
   C6  /   C5        C6  /   C5 
  /   /   / |       /  \.     | 
 /   /   /  |      /   .\     | 
*0  /  4R   |     /   /  4R   | 
|   |   |   |     |   |   |   | 
0   1   2   3     0   1   2   3 

 == Loci ==       == Locus 2 ==
 == 0 & 1==
```
From this we can infer that haplotype 0 has a mutation at position 0, haplotype 1 has a mutation at locus 2, haplotype 2 has a mutation at locus 1, and haplotype 3 has a mutation at loci 1 and 2.

The final line in the file (initial letter `S`) states that node 8 (the root) is inferred to have the sequence `000`, hence all mutations are from `0`→`1`, and the haplotypes are

```
        Haplotype
Locus |  0 1 2 3
------+---------
    0 |  1 0 0 0
    1 |  0 0 1 1
    2 |  0 1 0 1
```
which reconstructs the same data as in the original input file.
