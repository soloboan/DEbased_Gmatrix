DEbased_BlendedG <- function(DEresults,inpGeno,DEsample1,DEsample2,kbdist,wgtcolumn,s_tau,svalue_diagG,outname){
  library(MASS)
  cat('\n... import the DE Gene information ...\n')
  DE <- read.table(paste(DEresults,sep=''),header=T,stringsAsFactors=F)
  cat('... initial edits of DE data (e.g. exclude Log2FC = 0) ...\n')
  DE$CHR <- as.numeric(as.vector(substr(DE$contig_ID,4,6)))
  DE <- DE[order(DE$CHR,DE$gene_Start),]
  DE <- DE[which(DE$log2_FC!=0),]
  
  cat('... import the marker Position information ...\n')
  cat('    of the Genotype data - Plink map/bim file ...\n')
  mapinfo <- read.table(paste(inpGeno,'.bim',sep=''),stringsAsFactors = F)[,c(1,2,4)]
  colnames(mapinfo) <- c('CHR','MarkerName','Pos')
  pos <- mapinfo 
  cat('... initial edits (e.g. exclude zero chromosomes and unmmapped positions ...\n')
  pos <- pos[which(pos$Pos>0 & pos$CHR>0),]
  #table(paste(DE$sample1,DE$sample2,sep='+'))
  
  cat('\n... Using DE information of - sample ',DEsample1,' vs ',DEsample2,' ...\n')

  DE_ni <- DE[which(DE$sample1==DEsample1 & DE$sample2==DEsample2),]
  cat('...',nrow(DE_ni),' DE Genes to insert SNP markers ...\n')
  DE_markerinfo <- cbind.data.frame(pos[0,],DE_ni[0,])
  cat('... Positioning markers within ',kbdist/1000,' (kb) distance L/R of DE genes ...\n\n')
  
  iterchecks <- round(nrow(DE_ni)/5,digits=-1)
  for(strend in 1:nrow(DE_ni)){
    DE_snp <- pos[which(pos$CHR==DE_ni$CHR[strend] & 
                             pos$Pos>=DE_ni$gene_Start[strend]-kbdist & 
                             pos$Pos<=DE_ni$gene_END[strend]+kbdist),]
    if(nrow(DE_snp)>=1){
      DE_snp <- cbind.data.frame(DE_snp,DE_ni[strend,])
      DE_markerinfo <- rbind.data.frame(DE_markerinfo,DE_snp)  
    }
    if(strend %% iterchecks==0){
      cat('... ',strend,' ... out of ',nrow(DE_ni),' ...\n')
    }
  }
  #head(sort(table(DE_markerinfo$MarkerName),decreasing = T),25)
  #length(unique(DE_markerinfo$MarkerName))
  DEGene_ntimes <- data.frame(table(DE_markerinfo$MarkerName),stringsAsFactors=F)
  colnames(DEGene_ntimes) <- c('MarkerName','ntimes')
  DEGene_ntimesgt1 <- as.vector(DEGene_ntimes[which(DEGene_ntimes$ntimes>1),'MarkerName'])
  DEGene_markerinfo <- DE_markerinfo[!DE_markerinfo$MarkerName %in% DEGene_ntimesgt1,]
  cat('\n\n... ',nrow(DEGene_markerinfo),' Positioned in one (1) Gene ...\n')
  DEGene_ntimesgt1 <- DE_markerinfo[DE_markerinfo$MarkerName %in% DEGene_ntimesgt1,]
  
  DEG_gt1 <- unique(DEGene_ntimesgt1$MarkerName)
  cat('... ',length(DEG_gt1),' Positioned in multiple (>1) Gene ...\n')
  cat('\n... chosing one DE Gene for markers Positioned in several Genes ...\n')
  if(wgtcolumn=='log2_FC'){
    cat('... Choice based on maximum ',wgtcolumn,' value ...\n\n')  
  } else if (wgtcolumn=='p_value' | wgtcolumn=='q_value') {
    cat('... Choice based on minimum P-/Q- value ...\n\n')
  }
  iterchecks <- round(length(DEG_gt1)/5,digits=-1)
  for(DEMname in 1:length(DEG_gt1)){
    DE_snp <- DEGene_ntimesgt1[which(DEGene_ntimesgt1$MarkerName==DEG_gt1[DEMname]),]
    if(wgtcolumn=='log2_FC'){
      maxDE_snp <- DE_snp[which(abs(DE_snp[,wgtcolumn])==max(abs(DE_snp[,wgtcolumn]))),]  
      if(nrow(maxDE_snp)>1){
        maxDE_snp <- DE_snp[which(abs(DE_snp[,'q_value'])==min(abs(DE_snp[,'q_value']))),] 
      }
    } else if (wgtcolumn=='p_value' | wgtcolumn=='q_value') {
      maxDE_snp <- DE_snp[which(abs(DE_snp[,wgtcolumn])==min(abs(DE_snp[,wgtcolumn]))),]
      if(nrow(maxDE_snp)>1){
        maxDE_snp <- DE_snp[which(abs(DE_snp[,'log2_FC'])==max(abs(DE_snp[,'log2_FC']))),] 
      }
    }
    
    DEGene_markerinfo <- rbind.data.frame(DEGene_markerinfo,maxDE_snp)
    if(DEMname %% iterchecks==0){cat('... ',DEMname,' ... out of ',length(DEG_gt1),' ...\n')}
    rm(maxDE_snp)
  }
  cat('\n... ',nrow(DEGene_markerinfo),' (total) markers placed in DE Genes ...\n\n')
  cat('... exporting markers for building weighted G matrix ...\n')
  cat('     markers are stored in "DEGenes_markerinfo.txt" \n')
  write.table(DEGene_markerinfo$MarkerName,'DEGenes_markerinfo.txt',quote = F,row.names = F,col.names = F)
  
  cat('... Using PLINK extract markers for the DE and Polygenic G-matixes ...\n')
  mapinfo <- cbind.data.frame(mapinfo[,c('MarkerName')],AB=2)
  write.table(mapinfo,'recodeallele.txt',quote=F,row.names=F,col.names=F)
  system(paste('plink.exe --silent --chr-set 30 --nonfounders --bfile ',inpGeno,
               ' --extract DEGenes_markerinfo.txt --recode A --recode-allele recodeallele.txt',
               ' --out ',outname,'_DEGenes',sep=''))
  system(paste('plink.exe --silent --chr-set 30 --nonfounders --bfile ',inpGeno,
               ' --exclude DEGenes_markerinfo.txt --recode A --recode-allele recodeallele.txt',
               ' --out ',outname,'_polyg',sep=''))
  
  cat('\n\n... importing DE-based selected genotype markers ...\n')
  DEGenesMarkers <- read.table(paste(outname,'_DEGenes.raw',sep=''),stringsAsFactors=F,header=F,nrows=1)
  DEGenesMarkers <- t(DEGenesMarkers[,-1:-6])
  DEGenesMarkers <- gsub(pattern='_2',replacement='',x=DEGenesMarkers)
  colnames(DEGenesMarkers) <- c('MarkerName')
  DEGenesMarkers <- merge(DEGenesMarkers,DEGene_markerinfo,by='MarkerName',sort=F)
  DEGenesMarkers <- DEGenesMarkers[,c('MarkerName','log2_FC','p_value','q_value')]
  cat('... Building DE-based G-matrix ...\n')
  DEGenes <- read.table(paste(outname,'_DEGenes.raw',sep=''),stringsAsFactors=F,header=F,skip=1)
  idsgeno <- DEGenes[,1:2]
  DEGenes <- DEGenes[,-1:-6]
  cat('... computing allele frequency ...\n')
  freqA <- colMeans(DEGenes,na.rm = T)/2
  cat('... computing normalized weights ...\n')
  
  if(wgtcolumn=='log2_FC'){
    F2C <- abs(DEGenesMarkers[,wgtcolumn])
    wtnormalize <- F2C/mean(F2C)
  } else if(wgtcolumn=='p_value' | wgtcolumn=='q_value') {
    logPQ <- -log10(DEGenesMarkers[,wgtcolumn])
    wtnormalize <- logPQ/mean(logPQ)
  }
  
  #hist(wtnormalize,breaks=100)
  cat('... dealing with missing genotypes ...\n')
  DEGenes <- scale(DEGenes,center=T,scale=F)
  DEGenes[is.na(DEGenes)] <- 0
  cat('... Centering (-2p) and scaling matrix (using normalized wts)  ...\n')
  DEGenes <- scale(DEGenes,center=F,scale=wtnormalize)
  K <- round(sum(2*freqA*(1-freqA)),4)
  G_DE <- (tcrossprod(DEGenes))/K
  cat("... DE-based G-matrix computed (ZDZ'/K)  ...\n")
  
  cat('\n\n... importing Polygenic (non DE-based) genotype markers ...\n')
  poly <- read.table(paste(outname,'_polyg.raw',sep=''),stringsAsFactors = F,header = F,skip = 1)
  cat('... Building Polygenic (non DE-based) G-matrix ...\n')
  poly <- poly[,-1:-6]
  cat('... computing allele frequency ...\n')
  freqA_poly <- colMeans(poly,na.rm=T)/2
  cat('... dealing with missing genotypes ...\n')
  cat('... Centering (-2p) matrix ...\n')
  poly <- scale(poly,center=T,scale=F)
  poly[is.na(poly)] <- 0
  K_poly <- round(sum(2*freqA_poly*(1-freqA_poly)),4)
  G_poly <- (tcrossprod(poly))/K_poly
  cat("... Polygenic (non DE-based) G-matrix computed (ZZ'/K)  ...\n")
  
  cat('\n Blending DE-based and Polygenic (non DE-based) G-matrix ...\n')
  cat('... ',s_tau*100,'% of DE-based G-matrix and \n   ',100-(s_tau*100),'% of Polygenic (non DE-based) G-matrix ...\n')
  G_DEploy <- s_tau*G_DE + (1-s_tau)*G_poly
  
  cat('... adding a small value (',svalue_diagG,') to diagonals of Blended G-matrix ...\n')
  diag(G_DEploy) <- diag(G_DEploy)+svalue_diagG
  
  cat('... Inverting Blended G-matrix ...\n')
  G_DEployinv <- solve(G_DEploy)
  
  cat('\n... Exporting G and G inverse in asreml format ...\n')
  Glist <- as.data.frame(which(row(G_DEploy)>=col(G_DEploy),arr.ind=TRUE))
  Glist$G <- G_DEploy[lower.tri(G_DEploy,diag=T)]
  Glist <- Glist[,c(2,1,3)]
  Glist <- Glist[order(Glist[,2],Glist[,1]),]
  G_asreml <- Glist[,c(2,1,3)]
  write.table(G_asreml,paste(outname,'.grm',sep=''),quote=F,row.names=F,col.names=F)
  cat('... G-matrix done  ...\n')
  
  Glist <- as.data.frame(which(row(G_DEployinv)>=col(G_DEployinv),arr.ind=TRUE))
  Glist$G <- G_DEployinv[lower.tri(G_DEployinv,diag=T)]
  Glist <- Glist[,c(2,1,3)]
  Glist <- Glist[order(Glist[,2],Glist[,1]),]
  Ginv_asreml <- Glist[,c(2,1,3)]
  write.table(Ginv_asreml,paste(outname,'.giv',sep=''),quote=F,row.names=F,col.names=F)
  cat('... G inverse - matrix done  ...\n')
  
  cat ('\n')
  cat('... all Processes done ...\n')
  Gmat <- list(IDs=idsgeno,G_DE=G_DE,G_poly=G_poly,
               G=G_DEploy,Ginv=G_DEployinv,
               G_asreml=G_asreml,Ginv_asreml=Ginv_asreml)
  unlink(recursive=T,x=c('*.log','*.nosex','*.raw','recodeallele.txt'))
  return(Gmat)
}
