function class=knnclassification()
conn=database('cloud','root','krd123');
sqlquery='select TIR_GLDV_lh4x4_min,TIR_GLDV_sd16x16,TIR_GLDV_asm4x4_sd,TIR_GLDV_cs4x4_max,VIS_GLDV_con4x4_max,TIR_GLDV_cs4x4_mean,VIS_SADH_lh4x4_min,VIS_RF_sre,VIS_GLDV_ent16x16,TIR_GLDV_asm16x16,VIS_GLDV_ent4x4_min,VIS_GLDV_cs4x4_max,VIS_SF_maxpl,VIS_SF_minpl,TIR_SADH_cor4x4_min from feature';
setdbprefs('DataReturnFormat','numeric');
results=zeros(48,15);
y=zeros(48,1);
results = fetch(conn,sqlquery);
y=['A';'A';'A';'A';'B';'B';'B';'B';'C';'C';'C';'C';'D';'D';'D';'D';'E';'E';'E';'E';'F';'F';'F';'F';'G';'G';'G';'G';'H';'H';'H';'H';'I';'I';'I';'I';'J';'J';'J';'J';'K';'K';'K';'K';'L';'L';'L';'L'];
u=ClassificationKNN.fit(results,y);
xnew=[4.2281e-04,22.5551, 0.0183, 1.1903, 1.4023,0.4017,7.2764e-04, 0.1374, 0.5292, 0.0223,0, 2.5335,20,17,0.995];
class=predict(u,xnew);
end
