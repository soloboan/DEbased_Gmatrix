### Creating G matrix based on Differential Expression information  
#### SNP markers are weighted based on their location to a gene and the Log2 Fold change in gene expression of those genes   
##### script requires plin for extract markers and converting alleles to genotypes

                          s*G_DE + (1-s)*G_polygenic  

                                G_DE = ZDZ'/K  

                                  Z=M-2p  
                            D=1/normalized(weights)  
                                K=2*sum(p*q)  

#### Weights in D (see above) are based on the folowing  
        1. absolute log2 fold change from the DE analysis  
        2. -log10 Pvalue from the DE analysis  
        3. -log10 qvalue from the DE analysis  
    

This script is specific for our data, so columns in the DE results file reflect that.  
If you intend to use this then check columns or adapt columns to this format  

      Gene information column       ==> "gene_ID\\contig_ID\\gene_Start\\gene_END\\
      Sample infromation column     ==> "sample1\\sample2\\"
      Weigthing factor columns (DE) ==> "log2_FC\\p_value\\q_value"


#### The script require the following argument  
    DEresults         ==> The file containing the results from the DE analysis  
    inpGeno           ==> Input prefix of the PLINK format genotype file (plink binary format)  
    DEsample1         ==> Which DE comparison you want to use (eg. control/naive vs dead). first  
    DEsample2         ==> Which DE comparison you want to use (eg. control/naive vs dead). second  
    kbdist            ==> The left/Right of the SNP marker (reflect - distance of marker to gene)  
    wgtcolumn         ==> The weight to use (eg. Log_2FC or p-value or q-value)  
    s_tau             ==> blend the weighted matrix to polygenic matrix (s\*G_D + 1-s*Gpolygenic)  
    svalue_diagG      ==> small value to add to diagonals of G (stability and inversion problem)  
    outname           ==> output name of matrix and asreml format G matrix and inverse  


#### output 

- The script will give you the G matrix for the Blended DE G and Polygynic G  
- R object (This is a list) with the 
    - Orginal IDs
    - DE based G matrix
    - G matrix based on non-DE based markers
    - Blended G matrix 
    - Blended G inverse
    - Asreml format (G)
    - Asreml format (G inverse)
    
    



