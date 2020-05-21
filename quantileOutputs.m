function [percDiffB,percDiffP] = quantileOutputs(meanB,meanP,bB,bP,params)

biomQuantile=params(1);
prodQuantile=params(2);

predictedB=bB(1)+bB(2)*biomQuantile+bB(3)*prodQuantile+bB(4)*biomQuantile*prodQuantile;
predictedP=bP(1)+bP(2)*biomQuantile+bP(3)*prodQuantile+bP(4)*biomQuantile*prodQuantile;

percDiffB=(meanB-predictedB)*100/meanB;
percDiffP=(meanP-predictedP)*100/meanP;