source('DEbased_BlendedGmatrix.R')

## excuting example file with 
### naive vs AGD infection score > 3
### using the q values from the DE analysis of the Genes
### markers could be 10 kbs away from the Gene
### Blend 90% of the DE Gmatrix with 10% of the Polygenic part
#### add 0.01 to the diagonals of the matrix before inversion
DEGmat_qvalue_tau90_Kkb10 <- DEbased_BlendedG(DEresults='../AGDbolaksgeno/AGD_List.sig.genes.txt',
                                                  inpGeno='../AGDbolaksgeno/AGD_GS',
                                                  DEsample1='n',DEsample2='i2_3',
                                                  kbdist=10000,
                                                  wgtcolumn='q_value',
                                                  s_tau=0.90,
                                                  svalue_diagG=0.01,
                                                  outname='G_qvalue_tau90_Kkb10')

hist(diag(DEGmat_qvalue_tau90_Kkb10$G),breaks=25,col = sample(colours(),1),
     main='',xlab=' Diagonals of G')



## excuting example file with 
### naive vs AGD infection score > 3
### using the log2_FC from the DE analysis of the Genes
### markers could be 25 kbs away from the Gene
### Blend 90% of the DE Gmatrix with 10% of the Polygenic part
#### add 0.01 to the diagonals of the matrix before inversion
DEGmat_log2FC_tau90_Kkb25 <- DEbased_BlendedG(DEresults='../AGDbolaksgeno/AGD_List.sig.genes.txt',
                                              inpGeno='../AGDbolaksgeno/AGD_GS',
                                              DEsample1='n',DEsample2='i2_3',
                                              kbdist=25000,
                                              wgtcolumn='log2_FC',
                                              s_tau=0.90,
                                              svalue_diagG=0.01,
                                              outname='G_log2FC_tau90_Kkb25')

hist(diag(DEGmat_log2FC_tau90_Kkb25$G),breaks=25,col = sample(colours(),1),
     main='',xlab=' Diagonals of G')
