%function savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV_yrs,gainBLV_yrs,BLVw_yrs,gainBLVw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,BLV2_yrs,gainBLV2_yrs,BLV2w_yrs,gainBLV2w_yrs,BLV3_yrs,gainBLV3_yrs,BLV3w_yrs,gainBLV3w_yrs,BLV4_yrs,gainBLV4_yrs,BLV4w_yrs,gainBLV4w_yrs,BLV5_yrs,gainBLV5_yrs,BLV5w_yrs,gainBLV5w_yrs,BLV6_yrs,gainBLV6_yrs,BLV6w_yrs,gainBLV6w_yrs,Btrans,dBtrans,gainBtrans,r,a,z,K,flag,raR2,r_T,K_T,K_T_ratio,r_T_ratio,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,r2,a2,z2,K2,flag2,raR22,r_T2,K_T2,K_T_ratio2,r_T_ratio2,r3,a3,z3,K3,flag3,raR23,r_T3,K_T3,K_T_ratio3,r_T_ratio3,r4,a4,z4,K4,flag4,raR24,r_T4,K_T4,K_T_ratio4,r_T_ratio4,r5,a5,z5,K5,flag5,raR25,r_T5,K_T5,K_T_ratio5,r_T_ratio5,r6,a6,z6,K6,flag6,raR26,r_T6,K_T6,K_T_ratio6,r_T_ratio6,fitCode,P)
function savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,Btrans,dBtrans,gainBtrans,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,fitCode,P)

%save(FoodWebFile, 'Z_yrs','B_yrs','gainB_yrs','gainZ_yrs','v_yrs','TE_yrs','PB_yrs','TLik_yrs','TLi_yrs','TLk_yrs','TLall_yrs','Zw_yrs','Bw_yrs','gainBw_yrs','gainZw_yrs','vw_yrs','TEw_yrs','PBw_yrs','TLikw_yrs','TLiw_yrs','TLkw_yrs','TLallw_yrs','BLV_yrs','gainBLV_yrs','BLVw_yrs','gainBLVw_yrs','BLV1_yrs','gainBLV1_yrs','BLV1w_yrs','gainBLV1w_yrs','BLV2_yrs','gainBLV2_yrs','BLV2w_yrs','gainBLV2w_yrs','BLV3_yrs','gainBLV3_yrs','BLV3w_yrs','gainBLV3w_yrs','BLV4_yrs','gainBLV4_yrs','BLV4w_yrs','gainBLV4w_yrs','BLV5_yrs','gainBLV5_yrs','BLV5w_yrs','gainBLV5w_yrs','BLV6_yrs','gainBLV6_yrs','BLV6w_yrs','gainBLV6w_yrs','Btrans','dBtrans','gainBtrans','r','a','z','K','flag','raR2','r_T','K_T','K_T_ratio','r_T_ratio','r1','a1','z1','K1','flag1','raR21','r_T1','K_T1','K_T_ratio1','r_T_ratio1','r2','a2','z2','K2','flag2','raR22','r_T2','K_T2','K_T_ratio2','r_T_ratio2','r3','a3','z3','K3','flag3','raR23','r_T3','K_T3','K_T_ratio3','r_T_ratio3','r4','a4','z4','K4','flag4','raR24','r_T4','K_T4','K_T_ratio4','r_T_ratio4','r5','a5','z5','K5','flag5','raR25','r_T5','K_T5','K_T_ratio5','r_T_ratio5','r6','a6','z6','K6','flag6','raR26','r_T6','K_T6','K_T_ratio6','r_T_ratio6','fitCode','P');
save(FoodWebFile,'Z_yrs','B_yrs','gainB_yrs','gainZ_yrs','v_yrs','TE_yrs','PB_yrs','TLik_yrs','TLi_yrs','TLk_yrs','TLall_yrs','Zw_yrs','Bw_yrs','gainBw_yrs','gainZw_yrs','vw_yrs','TEw_yrs','PBw_yrs','TLikw_yrs','TLiw_yrs','TLkw_yrs','TLallw_yrs','BLV1_yrs','gainBLV1_yrs','BLV1w_yrs','gainBLV1w_yrs','Btrans','dBtrans','gainBtrans','r1','a1','z1','K1','flag1','raR21','r_T1','K_T1','K_T_ratio1','r_T_ratio1','fitCode','P');

end