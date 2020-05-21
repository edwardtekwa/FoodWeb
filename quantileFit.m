function SS = quantileFit(meanB,meanP,bB,bP,params)

biomQuantile=params(1);
prodQuantile=params(2);

predictedB=bB(1)+bB(2)*biomQuantile+bB(3)*prodQuantile+bB(4)*biomQuantile*prodQuantile;
predictedP=bP(1)+bP(2)*biomQuantile+bP(3)*prodQuantile+bP(4)*biomQuantile*prodQuantile;

percDiffB=abs(meanB-predictedB)/meanB;
percDiffP=abs(meanP-predictedP)/meanP;

%sumPercentDiff=percDiffB+percDiffP;
SS=((meanB-predictedB)^2)/meanB+((meanP-predictedP)^2)/meanP;