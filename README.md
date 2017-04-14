## Creating G matrix based on Differential Expression information  
### SNP markers are weighted in the based on their location to a gene and the Log2 Fold change in gene expression    

s*G_DE + (1-s)*G_polygenic  

G_DE = ZDZ'/K  

Z=M-2p  
D=1/normalized(weights)  
K=2*sum(p*q)  

weights =  
    1. absolute log2 fold change from the DE analysis  
    2. -log10 Pvalue from the DE analysis  
    3. -log10 qvalue from the DE analysis  
    

SNP is pecific for our data, so columns in the DE results file reflect that.  
If you intend to use this then check columns or adapt columns to this format  

"gene_ID\\contig_ID\\gene_Start\\gene_END\\sample1\\sample2\\test_Status\\FPKM_sample1\\	FPKM_sample2\\log2_FC\\p_value\\q_value"



#### The script require the following argument  
DEresults         ==> The file containing the results from the DE analysis  
inpGeno           ==> Input prefix of the PLINK format genotype file (plink binary format)  
DEsample1         ==> Which DE comparison you want to use (eg. control/naive vs dead). first  
DEsample2         ==> Which DE comparison you want to use (eg. control/naive vs dead). second  
kbdist            ==> The left/Right of the SNP marker (reflect - distance of marker to gene)  
wgtcolumn         ==> The weight to use (eg. Log_2FC or p-value or q-value)  
s_tau=0.90        ==> blend the weighted matrix to polygenic matrix (s*G_D + 1-s*Gpolygenic)  
svalue_diagG      ==> small value to add to diagonals of G (stability and inversion problem)  
outname           ==> output name of matrix and asreml format G matrix and inverse  


### output 

- The script will give you the G matrix for the Blended DE G and Polygynic G  
- R object (This is a list) with the 
    - Orginal IDs
    - DE based G matrix
    - G matrix based on non-DE based markers
    - Blended G matrix 
    - Blended G inverse
    
    



