C1 = (sum(pInc.dvHPC.Coor.dHPC(:,3)<0.001)/length(pInc.dvHPC.Coor.dHPC(:,3)))*100;
C2 = (sum(pDec.dvHPC.Coor.dHPC(:,3)<0.001)/length(pDec.dvHPC.Coor.dHPC(:,3)))*100;
R1 = 100 - C1 - C2;

C3 = (sum(pInc.dvHPC.Uncoor.dHPC(:,3)<0.001)/length(pInc.dvHPC.Uncoor.dHPC(:,3)))*100;
C4 = (sum(pDec.dvHPC.Uncoor.dHPC(:,3)<0.001)/length(pDec.dvHPC.Uncoor.dHPC(:,3)))*100;
R2 = 100 - C3 - C4;


C5 = (sum(pInc.dvHPC.Coor.vHPC(:,3)<0.001)/length(pInc.dvHPC.Coor.vHPC(:,3)))*100;
C6 = (sum(pDec.dvHPC.Coor.vHPC(:,3)<0.001)/length(pDec.dvHPC.Coor.vHPC(:,3)))*100;
R3 = 100 - C5 - C6;

C7= (sum(pInc.dvHPC.Uncoor.vHPC(:,3)<0.001)/length(pInc.dvHPC.Uncoor.vHPC(:,3)))*100;
C8 = (sum(pDec.dvHPC.Uncoor.vHPC(:,3)<0.001)/length(pDec.dvHPC.Uncoor.vHPC(:,3)))*100;
R4 = 100 - C7 - C8;

subplot(221),pie([C1 , C2, R1] , {'Increased' , 'Decreased' , 'Not Modulated'})
subplot(222),pie([C3 , C4, R2] , {'Increased' , 'Decreased' , 'Not Modulated'})
subplot(223),pie([C5 , C6, R3] , {'Increased' , 'Decreased' , 'Not Modulated'})
subplot(224),pie([C7 , C8, R4] , {'Increased' , 'Decreased' , 'Not Modulated'})