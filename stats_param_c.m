
BiomassSS=sum((totBiomass(1:12,1)-totBiomass(1:12,3)).^2)
ProdSS=sum((totProd(1:12,1)-totProd(1:12,3)).^2)
BiomassMeanDiff=100*mean((totBiomass(1:12,3)-totBiomass(1:12,1))./totBiomass(1:12,1))
ProdMeanDiff=100*mean((totProd(1:12,3)-totProd(1:12,1))./totProd(1:12,1))