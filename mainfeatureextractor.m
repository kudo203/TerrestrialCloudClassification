function a=mainfeatureextractor()
a=0;
h5files = dir('*.h5'); 
numfiles = length(h5files);
data_VIS = cell(1, numfiles);
data_TIR = cell(1,numfiles);
conn=database('mainfeatures','root','krd123');
for k = 1:numfiles
    data_VIS{k} = h5read(h5files(k).name,'/VHRR/Image Data/VIS1');
    data_TIR{k} = h5read(h5files(k).name,'/VHRR/Image Data/TIR1');
    file_name=h5files(k).name;
    d=',';
    c=')';
    i='''';
    %retrieve the values for VIS_SF features and transfer to the database
    [VIS_SF_maxpl,VIS_SF_minpl,VIS_SF_meanpl,VIS_SF_medianpl,VIS_SF_modepl,VIS_SF_rangepl,VIS_SF_sdpl]=VIS_SF(data_VIS{k});
    [TIR_SF_maxpl,TIR_SF_minpl,TIR_SF_meanpl,TIR_SF_medianpl,TIR_SF_modepl,TIR_SF_rangepl,TIR_SF_sdpl]=TIR_SF(data_TIR{k});
    [VIS_RF_sre,VIS_RF_lre,VIS_RF_gln,VIS_RF_rln,VIS_RF_rp]=VIS_RF(data_VIS{k});
    [TIR_RF_sre,TIR_RF_lre,TIR_RF_gln,TIR_RF_rln,TIR_RF_rp]=TIR_RF(data_TIR{k});
    
    [VIS_GLDV_mean16x16,VIS_GLDV_mean4x4_max,VIS_GLDV_mean4x4_min,VIS_GLDV_mean4x4_mean,VIS_GLDV_mean4x4_sd,VIS_GLDV_sd16x16,VIS_GLDV_sd4x4_max,VIS_GLDV_sd4x4_min,VIS_GLDV_sd4x4_mean,VIS_GLDV_sd4x4_sd,VIS_GLDV_asm16x16,VIS_GLDV_asm4x4_max,VIS_GLDV_asm4x4_min,VIS_GLDV_asm4x4_mean,VIS_GLDV_asm4x4_sd,VIS_GLDV_ent16x16,VIS_GLDV_ent4x4_max,VIS_GLDV_ent4x4_min,VIS_GLDV_ent4x4_mean,VIS_GLDV_ent4x4_sd,VIS_GLDV_lh16x16,VIS_GLDV_lh4x4_max,VIS_GLDV_lh4x4_min,VIS_GLDV_lh4x4_mean,VIS_GLDV_lh4x4_sd,VIS_GLDV_con16x16,VIS_GLDV_con4x4_max,VIS_GLDV_con4x4_min,VIS_GLDV_con4x4_mean,VIS_GLDV_con4x4_sd,VIS_GLDV_cs16x16,VIS_GLDV_cs4x4_max,VIS_GLDV_cs4x4_min,VIS_GLDV_cs4x4_mean,VIS_GLDV_cs4x4_sd,VIS_GLDV_cp16x16,VIS_GLDV_cp4x4_max,VIS_GLDV_cp4x4_min,VIS_GLDV_cp4x4_mean,VIS_GLDV_cp4x4_sd]=VIS_GLDV(data_VIS{k});
    [TIR_GLDV_mean16x16,TIR_GLDV_mean4x4_max,TIR_GLDV_mean4x4_min,TIR_GLDV_mean4x4_mean,TIR_GLDV_mean4x4_sd,TIR_GLDV_sd16x16,TIR_GLDV_sd4x4_max,TIR_GLDV_sd4x4_min,TIR_GLDV_sd4x4_mean,TIR_GLDV_sd4x4_sd,TIR_GLDV_asm16x16,TIR_GLDV_asm4x4_max,TIR_GLDV_asm4x4_min,TIR_GLDV_asm4x4_mean,TIR_GLDV_asm4x4_sd,TIR_GLDV_ent16x16,TIR_GLDV_ent4x4_max,TIR_GLDV_ent4x4_min,TIR_GLDV_ent4x4_mean,TIR_GLDV_ent4x4_sd,TIR_GLDV_lh16x16,TIR_GLDV_lh4x4_max,TIR_GLDV_lh4x4_min,TIR_GLDV_lh4x4_mean,TIR_GLDV_lh4x4_sd,TIR_GLDV_con16x16,TIR_GLDV_con4x4_max,TIR_GLDV_con4x4_min,TIR_GLDV_con4x4_mean,TIR_GLDV_con4x4_sd,TIR_GLDV_cs16x16,TIR_GLDV_cs4x4_max,TIR_GLDV_cs4x4_min,TIR_GLDV_cs4x4_mean,TIR_GLDV_cs4x4_sd,TIR_GLDV_cp16x16,TIR_GLDV_cp4x4_max,TIR_GLDV_cp4x4_min,TIR_GLDV_cp4x4_mean,TIR_GLDV_cp4x4_sd]=VIS_GLDV(data_TIR{k});
    
    [VIS_SADH_mean16x16,VIS_SADH_mean4x4_max,VIS_SADH_mean4x4_min,VIS_SADH_mean4x4_mean,VIS_SADH_mean4x4_sd,VIS_SADH_sd16x16,VIS_SADH_sd4x4_max,VIS_SADH_sd4x4_min,VIS_SADH_sd4x4_mean,VIS_SADH_sd4x4_sd,VIS_SADH_asm16x16,VIS_SADH_asm4x4_max,VIS_SADH_asm4x4_min,VIS_SADH_asm4x4_mean,VIS_SADH_asm4x4_sd,VIS_SADH_con16x16,VIS_SADH_con4x4_max,VIS_SADH_con4x4_min,VIS_SADH_con4x4_mean,VIS_SADH_con4x4_sd,VIS_SADH_cor16x16,VIS_SADH_cor4x4_max,VIS_SADH_cor4x4_min,VIS_SADH_cor4x4_mean,VIS_SADH_cor4x4_sd,VIS_SADH_ent16x16,VIS_SADH_ent4x4_max,VIS_SADH_ent4x4_min,VIS_SADH_ent4x4_mean,VIS_SADH_ent4x4_sd,VIS_SADH_lh16x16,VIS_SADH_lh4x4_max,VIS_SADH_lh4x4_min,VIS_SADH_lh4x4_mean,VIS_SADH_lh4x4_sd,VIS_SADH_cs16x16,VIS_SADH_cs4x4_max,VIS_SADH_cs4x4_min,VIS_SADH_cs4x4_mean,VIS_SADH_cs4x4_sd,VIS_SADH_cp16x16,VIS_SADH_cp4x4_max,VIS_SADH_cp4x4_min,VIS_SADH_cp4x4_mean,VIS_SADH_cp4x4_sd]=VIS_SADH(data_VIS{k});
    [TIR_SADH_mean16x16,TIR_SADH_mean4x4_max,TIR_SADH_mean4x4_min,TIR_SADH_mean4x4_mean,TIR_SADH_mean4x4_sd,TIR_SADH_sd16x16,TIR_SADH_sd4x4_max,TIR_SADH_sd4x4_min,TIR_SADH_sd4x4_mean,TIR_SADH_sd4x4_sd,TIR_SADH_asm16x16,TIR_SADH_asm4x4_max,TIR_SADH_asm4x4_min,TIR_SADH_asm4x4_mean,TIR_SADH_asm4x4_sd,TIR_SADH_con16x16,TIR_SADH_con4x4_max,TIR_SADH_con4x4_min,TIR_SADH_con4x4_mean,TIR_SADH_con4x4_sd,TIR_SADH_cor16x16,TIR_SADH_cor4x4_max,TIR_SADH_cor4x4_min,TIR_SADH_cor4x4_mean,TIR_SADH_cor4x4_sd,TIR_SADH_ent16x16,TIR_SADH_ent4x4_max,TIR_SADH_ent4x4_min,TIR_SADH_ent4x4_mean,TIR_SADH_ent4x4_sd,TIR_SADH_lh16x16,TIR_SADH_lh4x4_max,TIR_SADH_lh4x4_min,TIR_SADH_lh4x4_mean,TIR_SADH_lh4x4_sd,TIR_SADH_cs16x16,TIR_SADH_cs4x4_max,TIR_SADH_cs4x4_min,TIR_SADH_cs4x4_mean,TIR_SADH_cs4x4_sd,TIR_SADH_cp16x16,TIR_SADH_cp4x4_max,TIR_SADH_cp4x4_min,TIR_SADH_cp4x4_mean,TIR_SADH_cp4x4_sd]=TIR_SADH(data_TIR{k});
  
    [PF_IRCF,PF_LCF,PF_MCF,PF_CCF,PF_MCI,PF_CTT,PF_CA,PF_ST,PF_VICF] =PF();
    
    VIS_SF_str='insert into feature(file_name,VIS_SF_maxpl,VIS_SF_minpl,VIS_SF_meanpl,VIS_SF_medianpl,VIS_SF_modepl,VIS_SF_rangepl,VIS_SF_sdpl,';
    TIR_SF_str='TIR_SF_maxpl,TIR_SF_minpl,TIR_SF_meanpl,TIR_SF_medianpl,TIR_SF_modepl,TIR_SF_rangepl,TIR_SF_sdpl,';
    
    VIS_RF_str='VIS_RF_sre,VIS_RF_lre,VIS_RF_gln,VIS_RF_rp,VIS_RF_rln,';
    TIR_RF_str='TIR_RF_sre,TIR_RF_lre,TIR_RF_gln,TIR_RF_rp,TIR_RF_rln,';
    
    VIS_GLDV_str='VIS_GLDV_mean16x16,VIS_GLDV_mean4x4_max,VIS_GLDV_mean4x4_min,VIS_GLDV_mean4x4_mean,VIS_GLDV_mean4x4_sd,VIS_GLDV_sd16x16,VIS_GLDV_sd4x4_max,VIS_GLDV_sd4x4_min,VIS_GLDV_sd4x4_mean,VIS_GLDV_sd4x4_sd,VIS_GLDV_asm16x16,VIS_GLDV_asm4x4_max,VIS_GLDV_asm4x4_min,VIS_GLDV_asm4x4_mean,VIS_GLDV_asm4x4_sd,VIS_GLDV_ent16x16,VIS_GLDV_ent4x4_max,VIS_GLDV_ent4x4_min,VIS_GLDV_ent4x4_mean,VIS_GLDV_ent4x4_sd,VIS_GLDV_lh16x16,VIS_GLDV_lh4x4_max,VIS_GLDV_lh4x4_min,VIS_GLDV_lh4x4_mean,VIS_GLDV_lh4x4_sd,VIS_GLDV_con16x16,VIS_GLDV_con4x4_max,VIS_GLDV_con4x4_min,VIS_GLDV_con4x4_mean,VIS_GLDV_con4x4_sd,VIS_GLDV_cs16x16,VIS_GLDV_cs4x4_max,VIS_GLDV_cs4x4_min,VIS_GLDV_cs4x4_mean,VIS_GLDV_cs4x4_sd,VIS_GLDV_cp16x16,VIS_GLDV_cp4x4_max,VIS_GLDV_cp4x4_min,VIS_GLDV_cp4x4_mean,VIS_GLDV_cp4x4_sd,';
    TIR_GLDV_str='TIR_GLDV_mean16x16,TIR_GLDV_mean4x4_max,TIR_GLDV_mean4x4_min,TIR_GLDV_mean4x4_mean,TIR_GLDV_mean4x4_sd,TIR_GLDV_sd16x16,TIR_GLDV_sd4x4_max,TIR_GLDV_sd4x4_min,TIR_GLDV_sd4x4_mean,TIR_GLDV_sd4x4_sd,TIR_GLDV_asm16x16,TIR_GLDV_asm4x4_max,TIR_GLDV_asm4x4_min,TIR_GLDV_asm4x4_mean,TIR_GLDV_asm4x4_sd,TIR_GLDV_ent16x16,TIR_GLDV_ent4x4_max,TIR_GLDV_ent4x4_min,TIR_GLDV_ent4x4_mean,TIR_GLDV_ent4x4_sd,TIR_GLDV_lh16x16,TIR_GLDV_lh4x4_max,TIR_GLDV_lh4x4_min,TIR_GLDV_lh4x4_mean,TIR_GLDV_lh4x4_sd,TIR_GLDV_con16x16,TIR_GLDV_con4x4_max,TIR_GLDV_con4x4_min,TIR_GLDV_con4x4_mean,TIR_GLDV_con4x4_sd,TIR_GLDV_cs16x16,TIR_GLDV_cs4x4_max,TIR_GLDV_cs4x4_min,TIR_GLDV_cs4x4_mean,TIR_GLDV_cs4x4_sd,TIR_GLDV_cp16x16,TIR_GLDV_cp4x4_max,TIR_GLDV_cp4x4_min,TIR_GLDV_cp4x4_mean,TIR_GLDV_cp4x4_sd,';
    
    VIS_SADH_str='VIS_SADH_mean16x16,VIS_SADH_mean4x4_max,VIS_SADH_mean4x4_min,VIS_SADH_mean4x4_mean,VIS_SADH_mean4x4_sd,VIS_SADH_sd16x16,VIS_SADH_sd4x4_max,VIS_SADH_sd4x4_min,VIS_SADH_sd4x4_mean,VIS_SADH_sd4x4_sd,VIS_SADH_asm16x16,VIS_SADH_asm4x4_max,VIS_SADH_asm4x4_min,VIS_SADH_asm4x4_mean,VIS_SADH_asm4x4_sd,VIS_SADH_con16x16,VIS_SADH_con4x4_max,VIS_SADH_con4x4_min,VIS_SADH_con4x4_mean,VIS_SADH_con4x4_sd,VIS_SADH_cor16x16,VIS_SADH_cor4x4_max,VIS_SADH_cor4x4_min,VIS_SADH_cor4x4_mean,VIS_SADH_cor4x4_sd,VIS_SADH_ent16x16,VIS_SADH_ent4x4_max,VIS_SADH_ent4x4_min,VIS_SADH_ent4x4_mean,VIS_SADH_ent4x4_sd,VIS_SADH_lh16x16,VIS_SADH_lh4x4_max,VIS_SADH_lh4x4_min,VIS_SADH_lh4x4_mean,VIS_SADH_lh4x4_sd,VIS_SADH_cs16x16,VIS_SADH_cs4x4_max,VIS_SADH_cs4x4_min,VIS_SADH_cs4x4_mean,VIS_SADH_cs4x4_sd,VIS_SADH_cp16x16,VIS_SADH_cp4x4_max,VIS_SADH_cp4x4_min,VIS_SADH_cp4x4_mean,VIS_SADH_cp4x4_sd,';
    TIR_SADH_str='TIR_SADH_mean16x16,TIR_SADH_mean4x4_max,TIR_SADH_mean4x4_min,TIR_SADH_mean4x4_mean,TIR_SADH_mean4x4_sd,TIR_SADH_sd16x16,TIR_SADH_sd4x4_max,TIR_SADH_sd4x4_min,TIR_SADH_sd4x4_mean,TIR_SADH_sd4x4_sd,TIR_SADH_asm16x16,TIR_SADH_asm4x4_max,TIR_SADH_asm4x4_min,TIR_SADH_asm4x4_mean,TIR_SADH_asm4x4_sd,TIR_SADH_con16x16,TIR_SADH_con4x4_max,TIR_SADH_con4x4_min,TIR_SADH_con4x4_mean,TIR_SADH_con4x4_sd,TIR_SADH_cor16x16,TIR_SADH_cor4x4_max,TIR_SADH_cor4x4_min,TIR_SADH_cor4x4_mean,TIR_SADH_cor4x4_sd,TIR_SADH_ent16x16,TIR_SADH_ent4x4_max,TIR_SADH_ent4x4_min,TIR_SADH_ent4x4_mean,TIR_SADH_ent4x4_sd,TIR_SADH_lh16x16,TIR_SADH_lh4x4_max,TIR_SADH_lh4x4_min,TIR_SADH_lh4x4_mean,TIR_SADH_lh4x4_sd,TIR_SADH_cs16x16,TIR_SADH_cs4x4_max,TIR_SADH_cs4x4_min,TIR_SADH_cs4x4_mean,TIR_SADH_cs4x4_sd,TIR_SADH_cp16x16,TIR_SADH_cp4x4_max,TIR_SADH_cp4x4_min,TIR_SADH_cp4x4_mean,TIR_SADH_cp4x4_sd,';   
    PF_str ='PF_IRCF,PF_LCF,PF_MCF,PF_CCF,PF_MCI,PF_CTT,PF_CA,PF_ST,PF_VISCF)values(';
    VIS_SADH_cor4x4_max = rand();
    VIS_SADH_cor4x4_min = rand();
    VIS_SADH_cor4x4_mean = rand();
    VIS_SADH_cor4x4_sd = rand();
    TIR_SADH_cor4x4_max = rand();
    TIR_SADH_cor4x4_min = rand();
    TIR_SADH_cor4x4_mean = rand();
    TIR_SADH_cor4x4_sd = rand();
    
    VIS_SF_maxpl_val=num2str(VIS_SF_maxpl);
    VIS_SF_minpl_val=num2str(VIS_SF_minpl);
    VIS_SF_meanpl_val=num2str(VIS_SF_meanpl);
    VIS_SF_medianpl_val=num2str(VIS_SF_medianpl);
    VIS_SF_modepl_val=num2str(VIS_SF_modepl);
    VIS_SF_rangepl_val=num2str(VIS_SF_rangepl);
    VIS_SF_sdpl_val=num2str(VIS_SF_sdpl);
    
    TIR_SF_maxpl_val=num2str(TIR_SF_maxpl);
    TIR_SF_minpl_val=num2str(TIR_SF_minpl);
    TIR_SF_meanpl_val=num2str(TIR_SF_meanpl);
    TIR_SF_medianpl_val=num2str(TIR_SF_medianpl);
    TIR_SF_modepl_val=num2str(TIR_SF_modepl);
    TIR_SF_rangepl_val=num2str(TIR_SF_rangepl);
    TIR_SF_sdpl_val=num2str(TIR_SF_sdpl);
    
    VIS_RF_sre_val=num2str(VIS_RF_sre);
    VIS_RF_lre_val=num2str(VIS_RF_lre);
    VIS_RF_gln_val=num2str(VIS_RF_gln);
    VIS_RF_rln_val=num2str(VIS_RF_rln);
    VIS_RF_rp_val=num2str(VIS_RF_rp);
    
    TIR_RF_sre_val=num2str(TIR_RF_sre);
    TIR_RF_lre_val=num2str(TIR_RF_lre);
    TIR_RF_gln_val=num2str(TIR_RF_gln);
    TIR_RF_rln_val=num2str(TIR_RF_rln);
    TIR_RF_rp_val=num2str(TIR_RF_rp);
    
    VIS_GLDV_mean16x16_val=num2str(VIS_GLDV_mean16x16);
    VIS_GLDV_mean4x4_max_val=num2str(VIS_GLDV_mean4x4_max);
    VIS_GLDV_mean4x4_min_val=num2str(VIS_GLDV_mean4x4_min);
    VIS_GLDV_mean4x4_mean_val=num2str(VIS_GLDV_mean4x4_mean);
    VIS_GLDV_mean4x4_sd_val=num2str(VIS_GLDV_mean4x4_sd);
    VIS_GLDV_sd16x16_val=num2str(VIS_GLDV_sd16x16);
    VIS_GLDV_sd4x4_max_val=num2str(VIS_GLDV_sd4x4_max);
    VIS_GLDV_sd4x4_min_val=num2str(VIS_GLDV_sd4x4_min);
    VIS_GLDV_sd4x4_mean_val=num2str(VIS_GLDV_sd4x4_mean);
    VIS_GLDV_sd4x4_sd_val=num2str(VIS_GLDV_sd4x4_sd);
    VIS_GLDV_asm16x16_val=num2str(VIS_GLDV_asm16x16);
    VIS_GLDV_asm4x4_max_val=num2str(VIS_GLDV_asm4x4_max);
    VIS_GLDV_asm4x4_min_val=num2str(VIS_GLDV_asm4x4_min);
    VIS_GLDV_asm4x4_mean_val=num2str(VIS_GLDV_asm4x4_mean);
    VIS_GLDV_asm4x4_sd_val=num2str(VIS_GLDV_asm4x4_sd);
    VIS_GLDV_ent16x16_val=num2str(VIS_GLDV_ent16x16);
    VIS_GLDV_ent4x4_max_val=num2str(VIS_GLDV_ent4x4_max);
    VIS_GLDV_ent4x4_min_val=num2str(VIS_GLDV_ent4x4_min);
    VIS_GLDV_ent4x4_mean_val=num2str(VIS_GLDV_ent4x4_mean);
    VIS_GLDV_ent4x4_sd_val=num2str(VIS_GLDV_ent4x4_sd);
    VIS_GLDV_lh16x16_val=num2str(VIS_GLDV_lh16x16);
    VIS_GLDV_lh4x4_max_val=num2str(VIS_GLDV_lh4x4_max);
    VIS_GLDV_lh4x4_min_val=num2str(VIS_GLDV_lh4x4_min);
    VIS_GLDV_lh4x4_mean_val=num2str(VIS_GLDV_lh4x4_mean);
    VIS_GLDV_lh4x4_sd_val=num2str(VIS_GLDV_lh4x4_sd);
    VIS_GLDV_con16x16_val=num2str(VIS_GLDV_con16x16);
    VIS_GLDV_con4x4_max_val=num2str(VIS_GLDV_con4x4_max);
    VIS_GLDV_con4x4_min_val=num2str(VIS_GLDV_con4x4_min);
    VIS_GLDV_con4x4_mean_val=num2str(VIS_GLDV_con4x4_mean);
    VIS_GLDV_con4x4_sd_val=num2str(VIS_GLDV_con4x4_sd);
    VIS_GLDV_cs16x16_val=num2str(VIS_GLDV_cs16x16);
    VIS_GLDV_cs4x4_max_val=num2str(VIS_GLDV_cs4x4_max);
    VIS_GLDV_cs4x4_min_val=num2str(VIS_GLDV_cs4x4_min);
    VIS_GLDV_cs4x4_mean_val=num2str(VIS_GLDV_cs4x4_mean);
    VIS_GLDV_cs4x4_sd_val=num2str(VIS_GLDV_cs4x4_sd);
    VIS_GLDV_cp16x16_val=num2str(VIS_GLDV_cp16x16);
    VIS_GLDV_cp4x4_max_val=num2str(VIS_GLDV_cp4x4_max);
    VIS_GLDV_cp4x4_min_val=num2str(VIS_GLDV_cp4x4_min);
    VIS_GLDV_cp4x4_mean_val=num2str(VIS_GLDV_cp4x4_mean);
    VIS_GLDV_cp4x4_sd_val=num2str(VIS_GLDV_cp4x4_sd);
    
    TIR_GLDV_mean16x16_val=num2str(TIR_GLDV_mean16x16);
    TIR_GLDV_mean4x4_max_val=num2str(TIR_GLDV_mean4x4_max);
    TIR_GLDV_mean4x4_min_val=num2str(TIR_GLDV_mean4x4_min);
    TIR_GLDV_mean4x4_mean_val=num2str(TIR_GLDV_mean4x4_mean);
    TIR_GLDV_mean4x4_sd_val=num2str(TIR_GLDV_mean4x4_sd);
    TIR_GLDV_sd16x16_val=num2str(TIR_GLDV_sd16x16);
    TIR_GLDV_sd4x4_max_val=num2str(TIR_GLDV_sd4x4_max);
    TIR_GLDV_sd4x4_min_val=num2str(TIR_GLDV_sd4x4_min);
    TIR_GLDV_sd4x4_mean_val=num2str(TIR_GLDV_sd4x4_mean);
    TIR_GLDV_sd4x4_sd_val=num2str(TIR_GLDV_sd4x4_sd);
    TIR_GLDV_asm16x16_val=num2str(TIR_GLDV_asm16x16);
    TIR_GLDV_asm4x4_max_val=num2str(TIR_GLDV_asm4x4_max);
    TIR_GLDV_asm4x4_min_val=num2str(TIR_GLDV_asm4x4_min);
    TIR_GLDV_asm4x4_mean_val=num2str(TIR_GLDV_asm4x4_mean);
    TIR_GLDV_asm4x4_sd_val=num2str(TIR_GLDV_asm4x4_sd);
    TIR_GLDV_ent16x16_val=num2str(TIR_GLDV_ent16x16);
    TIR_GLDV_ent4x4_max_val=num2str(TIR_GLDV_ent4x4_max);
    TIR_GLDV_ent4x4_min_val=num2str(TIR_GLDV_ent4x4_min);
    TIR_GLDV_ent4x4_mean_val=num2str(TIR_GLDV_ent4x4_mean);
    TIR_GLDV_ent4x4_sd_val=num2str(TIR_GLDV_ent4x4_sd);
    TIR_GLDV_lh16x16_val=num2str(TIR_GLDV_lh16x16);
    TIR_GLDV_lh4x4_max_val=num2str(TIR_GLDV_lh4x4_max);
    TIR_GLDV_lh4x4_min_val=num2str(TIR_GLDV_lh4x4_min);
    TIR_GLDV_lh4x4_mean_val=num2str(TIR_GLDV_lh4x4_mean);
    TIR_GLDV_lh4x4_sd_val=num2str(TIR_GLDV_lh4x4_sd);
    TIR_GLDV_con16x16_val=num2str(TIR_GLDV_con16x16);
    TIR_GLDV_con4x4_max_val=num2str(TIR_GLDV_con4x4_max);
    TIR_GLDV_con4x4_min_val=num2str(TIR_GLDV_con4x4_min);
    TIR_GLDV_con4x4_mean_val=num2str(TIR_GLDV_con4x4_mean);
    TIR_GLDV_con4x4_sd_val=num2str(TIR_GLDV_con4x4_sd);
    TIR_GLDV_cs16x16_val=num2str(TIR_GLDV_cs16x16);
    TIR_GLDV_cs4x4_max_val=num2str(TIR_GLDV_cs4x4_max);
    TIR_GLDV_cs4x4_min_val=num2str(TIR_GLDV_cs4x4_min);
    TIR_GLDV_cs4x4_mean_val=num2str(TIR_GLDV_cs4x4_mean);
    TIR_GLDV_cs4x4_sd_val=num2str(TIR_GLDV_cs4x4_sd);
    TIR_GLDV_cp16x16_val=num2str(TIR_GLDV_cp16x16);
    TIR_GLDV_cp4x4_max_val=num2str(TIR_GLDV_cp4x4_max);
    TIR_GLDV_cp4x4_min_val=num2str(TIR_GLDV_cp4x4_min);
    TIR_GLDV_cp4x4_mean_val=num2str(TIR_GLDV_cp4x4_mean);
    TIR_GLDV_cp4x4_sd_val=num2str(TIR_GLDV_cp4x4_sd);
    
    VIS_SADH_mean16x16_val=num2str(VIS_SADH_mean16x16);
    VIS_SADH_mean4x4_max_val=num2str(VIS_SADH_mean4x4_max);
    VIS_SADH_mean4x4_min_val=num2str(VIS_SADH_mean4x4_min);
    VIS_SADH_mean4x4_mean_val=num2str(VIS_SADH_mean4x4_mean);
    VIS_SADH_mean4x4_sd_val=num2str(VIS_SADH_mean4x4_sd);
    VIS_SADH_sd16x16_val=num2str(VIS_SADH_sd16x16);
    VIS_SADH_sd4x4_max_val=num2str(VIS_SADH_sd4x4_max);
    VIS_SADH_sd4x4_min_val=num2str(VIS_SADH_sd4x4_min);
    VIS_SADH_sd4x4_mean_val=num2str(VIS_SADH_sd4x4_mean);
    VIS_SADH_sd4x4_sd_val=num2str(VIS_SADH_sd4x4_sd);
    VIS_SADH_asm16x16_val=num2str(VIS_SADH_asm16x16);
    VIS_SADH_asm4x4_max_val=num2str(VIS_SADH_asm4x4_max);
    VIS_SADH_asm4x4_min_val=num2str(VIS_SADH_asm4x4_min);
    VIS_SADH_asm4x4_mean_val=num2str(VIS_SADH_asm4x4_mean);
    VIS_SADH_asm4x4_sd_val=num2str(VIS_SADH_asm4x4_sd);
    VIS_SADH_ent16x16_val=num2str(VIS_SADH_ent16x16);
    VIS_SADH_ent4x4_max_val=num2str(VIS_SADH_ent4x4_max);
    VIS_SADH_ent4x4_min_val=num2str(VIS_SADH_ent4x4_min);
    VIS_SADH_ent4x4_mean_val=num2str(VIS_SADH_ent4x4_mean);
    VIS_SADH_ent4x4_sd_val=num2str(VIS_SADH_ent4x4_sd);
    VIS_SADH_lh16x16_val=num2str(VIS_SADH_lh16x16);
    VIS_SADH_lh4x4_max_val=num2str(VIS_SADH_lh4x4_max);
    VIS_SADH_lh4x4_min_val=num2str(VIS_SADH_lh4x4_min);
    VIS_SADH_lh4x4_mean_val=num2str(VIS_SADH_lh4x4_mean);
    VIS_SADH_lh4x4_sd_val=num2str(VIS_SADH_lh4x4_sd);
    VIS_SADH_con16x16_val=num2str(VIS_SADH_con16x16);
    VIS_SADH_con4x4_max_val=num2str(VIS_SADH_con4x4_max);
    VIS_SADH_con4x4_min_val=num2str(VIS_SADH_con4x4_min);
    VIS_SADH_con4x4_mean_val=num2str(VIS_SADH_con4x4_mean);
    VIS_SADH_con4x4_sd_val=num2str(VIS_SADH_con4x4_sd);
    VIS_SADH_cor16x16_val=num2str(VIS_SADH_cor16x16);
    VIS_SADH_cor4x4_max_val=num2str(VIS_SADH_cor4x4_max);
    VIS_SADH_cor4x4_min_val=num2str(VIS_SADH_cor4x4_min);
    VIS_SADH_cor4x4_mean_val=num2str(VIS_SADH_cor4x4_mean);
    VIS_SADH_cor4x4_sd_val=num2str(VIS_SADH_cor4x4_sd);
    VIS_SADH_cs16x16_val=num2str(VIS_SADH_cs16x16);
    VIS_SADH_cs4x4_max_val=num2str(VIS_SADH_cs4x4_max);
    VIS_SADH_cs4x4_min_val=num2str(VIS_SADH_cs4x4_min);
    VIS_SADH_cs4x4_mean_val=num2str(VIS_SADH_cs4x4_mean);
    VIS_SADH_cs4x4_sd_val=num2str(VIS_SADH_cs4x4_sd);
    VIS_SADH_cp16x16_val=num2str(VIS_SADH_cp16x16);
    VIS_SADH_cp4x4_max_val=num2str(VIS_SADH_cp4x4_max);
    VIS_SADH_cp4x4_min_val=num2str(VIS_SADH_cp4x4_min);
    VIS_SADH_cp4x4_mean_val=num2str(VIS_SADH_cp4x4_mean);
    VIS_SADH_cp4x4_sd_val=num2str(VIS_SADH_cp4x4_sd);
    
    TIR_SADH_mean16x16_val=num2str(TIR_SADH_mean16x16);
    TIR_SADH_mean4x4_max_val=num2str(TIR_SADH_mean4x4_max);
    TIR_SADH_mean4x4_min_val=num2str(TIR_SADH_mean4x4_min);
    TIR_SADH_mean4x4_mean_val=num2str(TIR_SADH_mean4x4_mean);
    TIR_SADH_mean4x4_sd_val=num2str(TIR_SADH_mean4x4_sd);
    TIR_SADH_sd16x16_val=num2str(TIR_SADH_sd16x16);
    TIR_SADH_sd4x4_max_val=num2str(TIR_SADH_sd4x4_max);
    TIR_SADH_sd4x4_min_val=num2str(TIR_SADH_sd4x4_min);
    TIR_SADH_sd4x4_mean_val=num2str(TIR_SADH_sd4x4_mean);
    TIR_SADH_sd4x4_sd_val=num2str(TIR_SADH_sd4x4_sd);
    TIR_SADH_asm16x16_val=num2str(TIR_SADH_asm16x16);
    TIR_SADH_asm4x4_max_val=num2str(TIR_SADH_asm4x4_max);
    TIR_SADH_asm4x4_min_val=num2str(TIR_SADH_asm4x4_min);
    TIR_SADH_asm4x4_mean_val=num2str(TIR_SADH_asm4x4_mean);
    TIR_SADH_asm4x4_sd_val=num2str(TIR_SADH_asm4x4_sd);
    TIR_SADH_ent16x16_val=num2str(TIR_SADH_ent16x16);
    TIR_SADH_ent4x4_max_val=num2str(TIR_SADH_ent4x4_max);
    TIR_SADH_ent4x4_min_val=num2str(TIR_SADH_ent4x4_min);
    TIR_SADH_ent4x4_mean_val=num2str(TIR_SADH_ent4x4_mean);
    TIR_SADH_ent4x4_sd_val=num2str(TIR_SADH_ent4x4_sd);
    TIR_SADH_lh16x16_val=num2str(TIR_SADH_lh16x16);
    TIR_SADH_lh4x4_max_val=num2str(TIR_SADH_lh4x4_max);
    TIR_SADH_lh4x4_min_val=num2str(TIR_SADH_lh4x4_min);
    TIR_SADH_lh4x4_mean_val=num2str(TIR_SADH_lh4x4_mean);
    TIR_SADH_lh4x4_sd_val=num2str(TIR_SADH_lh4x4_sd);
    TIR_SADH_con16x16_val=num2str(TIR_SADH_con16x16);
    TIR_SADH_con4x4_max_val=num2str(TIR_SADH_con4x4_max);
    TIR_SADH_con4x4_min_val=num2str(TIR_SADH_con4x4_min);
    TIR_SADH_con4x4_mean_val=num2str(TIR_SADH_con4x4_mean);
    TIR_SADH_con4x4_sd_val=num2str(TIR_SADH_con4x4_sd);
    TIR_SADH_cor16x16_val=num2str(TIR_SADH_cor16x16);
    TIR_SADH_cor4x4_max_val=num2str(TIR_SADH_cor4x4_max);
    TIR_SADH_cor4x4_min_val=num2str(TIR_SADH_cor4x4_min);
    TIR_SADH_cor4x4_mean_val=num2str(TIR_SADH_cor4x4_mean);
    TIR_SADH_cor4x4_sd_val=num2str(TIR_SADH_cor4x4_sd);
    TIR_SADH_cs16x16_val=num2str(TIR_SADH_cs16x16);
    TIR_SADH_cs4x4_max_val=num2str(TIR_SADH_cs4x4_max);
    TIR_SADH_cs4x4_min_val=num2str(TIR_SADH_cs4x4_min);
    TIR_SADH_cs4x4_mean_val=num2str(TIR_SADH_cs4x4_mean);
    TIR_SADH_cs4x4_sd_val=num2str(TIR_SADH_cs4x4_sd);
    TIR_SADH_cp16x16_val=num2str(TIR_SADH_cp16x16);
    TIR_SADH_cp4x4_max_val=num2str(TIR_SADH_cp4x4_max);
    TIR_SADH_cp4x4_min_val=num2str(TIR_SADH_cp4x4_min);
    TIR_SADH_cp4x4_mean_val=num2str(TIR_SADH_cp4x4_mean);
    TIR_SADH_cp4x4_sd_val=num2str(TIR_SADH_cp4x4_sd);
    
    PF_IRCF_val=num2str(PF_IRCF);
    PF_LCF_val=num2str(PF_LCF);
    PF_MCF_val=num2str(PF_MCF);
    PF_CCF_val=num2str(PF_CCF);
    PF_MCI_val=num2str(PF_MCI);
    PF_CTT_val=num2str(PF_CTT);
    PF_CA_val=num2str(PF_CA);
    PF_ST_val=num2str(PF_ST);
    PF_VISCF_val=num2str(PF_VICF);
    
    initial_query='insert into features(TIR_GLDV_lh4x4_min)values(';
    value_query_VIS_SF=strcat(VIS_SF_maxpl_val,d,VIS_SF_minpl_val,d,VIS_SF_meanpl_val,d,VIS_SF_medianpl_val,d,VIS_SF_modepl_val,d,VIS_SF_rangepl_val,d,VIS_SF_sdpl_val);
    value_query_TIR_SF=strcat(TIR_SF_maxpl_val,d,TIR_SF_minpl_val,d,TIR_SF_meanpl_val,d,TIR_SF_medianpl_val,d,TIR_SF_modepl_val,d,TIR_SF_rangepl_val,d,TIR_SF_sdpl_val);
    
    value_query_VIS_RF=strcat(VIS_RF_sre_val,d,VIS_RF_lre_val,d,VIS_RF_gln_val,d,VIS_RF_rp_val,d,VIS_RF_rln_val);
    value_query_TIR_RF=strcat(TIR_RF_sre_val,d,TIR_RF_lre_val,d,TIR_RF_gln_val,d,TIR_RF_rp_val,d,TIR_RF_rln_val);
    
    value_query_VIS_GLDV=strcat(VIS_GLDV_mean16x16_val,d,VIS_GLDV_mean4x4_max_val,d,VIS_GLDV_mean4x4_min_val,d,VIS_GLDV_mean4x4_mean_val,d,VIS_GLDV_mean4x4_sd_val,d,VIS_GLDV_sd16x16_val,d,VIS_GLDV_sd4x4_max_val,d,VIS_GLDV_sd4x4_min_val,d,VIS_GLDV_sd4x4_mean_val,d,VIS_GLDV_sd4x4_sd_val,d,VIS_GLDV_asm16x16_val,d,VIS_GLDV_asm4x4_max_val,d,VIS_GLDV_asm4x4_min_val,d,VIS_GLDV_asm4x4_mean_val,d,VIS_GLDV_asm4x4_sd_val,d,VIS_GLDV_ent16x16_val,d,VIS_GLDV_ent4x4_max_val,d,VIS_GLDV_ent4x4_min_val,d,VIS_GLDV_ent4x4_mean_val,d,VIS_GLDV_ent4x4_sd_val,d,VIS_GLDV_lh16x16_val,d,VIS_GLDV_lh4x4_max_val,d,VIS_GLDV_lh4x4_min_val,d,VIS_GLDV_lh4x4_mean_val,d,VIS_GLDV_lh4x4_sd_val,d,VIS_GLDV_con16x16_val,d,VIS_GLDV_con4x4_max_val,d,VIS_GLDV_con4x4_min_val,d,VIS_GLDV_con4x4_mean_val,d,VIS_GLDV_con4x4_sd_val,d,VIS_GLDV_cs16x16_val,d,VIS_GLDV_cs4x4_max_val,d,VIS_GLDV_cs4x4_min_val,d,VIS_GLDV_cs4x4_mean_val,d,VIS_GLDV_cs4x4_sd_val,d,VIS_GLDV_cp16x16_val,d,VIS_GLDV_cp4x4_max_val,d,VIS_GLDV_cp4x4_min_val,d,VIS_GLDV_cp4x4_mean_val,d,VIS_GLDV_cp4x4_sd_val);
    value_query_TIR_GLDV=strcat(TIR_GLDV_mean16x16_val,d,TIR_GLDV_mean4x4_max_val,d,TIR_GLDV_mean4x4_min_val,d,TIR_GLDV_mean4x4_mean_val,d,TIR_GLDV_mean4x4_sd_val,d,TIR_GLDV_sd16x16_val,d,TIR_GLDV_sd4x4_max_val,d,TIR_GLDV_sd4x4_min_val,d,TIR_GLDV_sd4x4_mean_val,d,TIR_GLDV_sd4x4_sd_val,d,TIR_GLDV_asm16x16_val,d,TIR_GLDV_asm4x4_max_val,d,TIR_GLDV_asm4x4_min_val,d,TIR_GLDV_asm4x4_mean_val,d,TIR_GLDV_asm4x4_sd_val,d,TIR_GLDV_ent16x16_val,d,TIR_GLDV_ent4x4_max_val,d,TIR_GLDV_ent4x4_min_val,d,TIR_GLDV_ent4x4_mean_val,d,TIR_GLDV_ent4x4_sd_val,d,TIR_GLDV_lh16x16_val,d,TIR_GLDV_lh4x4_max_val,d,TIR_GLDV_lh4x4_min_val,d,TIR_GLDV_lh4x4_mean_val,d,TIR_GLDV_lh4x4_sd_val,d,TIR_GLDV_con16x16_val,d,TIR_GLDV_con4x4_max_val,d,TIR_GLDV_con4x4_min_val,d,TIR_GLDV_con4x4_mean_val,d,TIR_GLDV_con4x4_sd_val,d,TIR_GLDV_cs16x16_val,d,TIR_GLDV_cs4x4_max_val,d,TIR_GLDV_cs4x4_min_val,d,TIR_GLDV_cs4x4_mean_val,d,TIR_GLDV_cs4x4_sd_val,d,TIR_GLDV_cp16x16_val,d,TIR_GLDV_cp4x4_max_val,d,TIR_GLDV_cp4x4_min_val,d,TIR_GLDV_cp4x4_mean_val,d,TIR_GLDV_cp4x4_sd_val);
    
    value_query_VIS_SADH=strcat(VIS_SADH_mean16x16_val,d,VIS_SADH_mean4x4_max_val,d,VIS_SADH_mean4x4_min_val,d,VIS_SADH_mean4x4_mean_val,d,VIS_SADH_mean4x4_sd_val,d,VIS_SADH_sd16x16_val,d,VIS_SADH_sd4x4_max_val,d,VIS_SADH_sd4x4_min_val,d,VIS_SADH_sd4x4_mean_val,d,VIS_SADH_sd4x4_sd_val,d,VIS_SADH_asm16x16_val,d,VIS_SADH_asm4x4_max_val,d,VIS_SADH_asm4x4_min_val,d,VIS_SADH_asm4x4_mean_val,d,VIS_SADH_asm4x4_sd_val,d,VIS_SADH_con16x16_val,d,VIS_SADH_con4x4_max_val,d,VIS_SADH_con4x4_min_val,d,VIS_SADH_con4x4_mean_val,d,VIS_SADH_con4x4_sd_val,d,VIS_SADH_cor16x16_val,d,VIS_SADH_cor4x4_max_val,d,VIS_SADH_cor4x4_min_val,d,VIS_SADH_cor4x4_mean_val,d,VIS_SADH_cor4x4_sd_val,d,VIS_SADH_ent16x16_val,d,VIS_SADH_ent4x4_max_val,d,VIS_SADH_ent4x4_min_val,d,VIS_SADH_ent4x4_mean_val,d,VIS_SADH_ent4x4_sd_val,d,VIS_SADH_lh16x16_val,d,VIS_SADH_lh4x4_max_val,d,VIS_SADH_lh4x4_min_val,d,VIS_SADH_lh4x4_mean_val,d,VIS_SADH_lh4x4_sd_val,d,VIS_SADH_cs16x16_val,d,VIS_SADH_cs4x4_max_val,d,VIS_SADH_cs4x4_min_val,d,VIS_SADH_cs4x4_mean_val,d,VIS_SADH_cs4x4_sd_val,d,VIS_SADH_cp16x16_val,d,VIS_SADH_cp4x4_max_val,d,VIS_SADH_cp4x4_min_val,d,VIS_SADH_cp4x4_mean_val,d,VIS_SADH_cp4x4_sd_val);
    value_query_TIR_SADH1=strcat(TIR_SADH_mean16x16_val,d,TIR_SADH_mean4x4_max_val,d,TIR_SADH_mean4x4_min_val,d,TIR_SADH_mean4x4_mean_val,d,TIR_SADH_mean4x4_sd_val,d,TIR_SADH_sd16x16_val,d,TIR_SADH_sd4x4_max_val,d,TIR_SADH_sd4x4_min_val,d,TIR_SADH_sd4x4_mean_val,d,TIR_SADH_sd4x4_sd_val,d,TIR_SADH_asm16x16_val,d,TIR_SADH_asm4x4_max_val,d,TIR_SADH_asm4x4_min_val,d,TIR_SADH_asm4x4_mean_val,d,TIR_SADH_asm4x4_sd_val,d,TIR_SADH_con16x16_val,d,TIR_SADH_con4x4_max_val,d,TIR_SADH_con4x4_min_val,d,TIR_SADH_con4x4_mean_val,d,TIR_SADH_con4x4_sd_val,d,TIR_SADH_cor16x16_val,d,TIR_SADH_cor4x4_max_val,d,TIR_SADH_cor4x4_min_val,d,TIR_SADH_cor4x4_mean_val,d,TIR_SADH_cor4x4_sd_val,d,TIR_SADH_ent16x16_val,d,TIR_SADH_ent4x4_max_val,d,TIR_SADH_ent4x4_min_val,d,TIR_SADH_ent4x4_mean_val,d,TIR_SADH_ent4x4_sd_val,d,TIR_SADH_lh16x16_val,d,TIR_SADH_lh4x4_max_val,d,TIR_SADH_lh4x4_min_val,d,TIR_SADH_lh4x4_mean_val,d,TIR_SADH_lh4x4_sd_val,d,TIR_SADH_cs16x16_val,d,TIR_SADH_cs4x4_max_val,d,TIR_SADH_cs4x4_min_val,d,TIR_SADH_cs4x4_mean_val,d,TIR_SADH_cs4x4_sd_val,d,TIR_SADH_cp16x16_val,d,TIR_SADH_cp4x4_max_val,d,TIR_SADH_cp4x4_min_val,d,TIR_SADH_cp4x4_mean_val,d,TIR_SADH_cp4x4_sd_val);
    
    value_query_TIR_SADH=strcat(TIR_SADH_mean16x16_val,d,TIR_SADH_mean4x4_max_val,d,TIR_SADH_mean4x4_min_val,d,TIR_SADH_mean4x4_mean_val,d,TIR_SADH_mean4x4_sd_val,d,TIR_SADH_sd16x16_val,d,TIR_SADH_sd4x4_max_val,d,TIR_SADH_sd4x4_min_val,d,TIR_SADH_sd4x4_mean_val,d,TIR_SADH_sd4x4_sd_val,d,TIR_SADH_asm16x16_val,d,TIR_SADH_asm4x4_max_val,d,TIR_SADH_asm4x4_min_val,d,TIR_SADH_asm4x4_mean_val,d,TIR_SADH_asm4x4_sd_val,d,TIR_SADH_con16x16_val,d,TIR_SADH_con4x4_max_val,d,TIR_SADH_con4x4_min_val,d,TIR_SADH_con4x4_mean_val,d,TIR_SADH_con4x4_sd_val,d,TIR_SADH_cor16x16_val,d,TIR_SADH_ent16x16_val,d,TIR_SADH_ent4x4_max_val,d,TIR_SADH_ent4x4_min_val,d,TIR_SADH_ent4x4_mean_val,d,TIR_SADH_ent4x4_sd_val,d,TIR_SADH_lh16x16_val,d,TIR_SADH_lh4x4_max_val,d,TIR_SADH_lh4x4_min_val,d,TIR_SADH_lh4x4_mean_val,d,TIR_SADH_lh4x4_sd_val,d,TIR_SADH_cs16x16_val,d,TIR_SADH_cs4x4_max_val,d,TIR_SADH_cs4x4_min_val,d,TIR_SADH_cs4x4_mean_val,d,TIR_SADH_cs4x4_sd_val,d,TIR_SADH_cp16x16_val,d,TIR_SADH_cp4x4_max_val,d,TIR_SADH_cp4x4_min_val,d,TIR_SADH_cp4x4_mean_val,d,TIR_SADH_cp4x4_sd_val);
    
    value_query_PF=strcat(PF_IRCF_val,d,PF_LCF_val,d,PF_MCF_val,d,PF_CCF_val,d,PF_MCI_val,d,PF_CTT_val,d,PF_CA_val,d,PF_ST_val,d,PF_VISCF_val);
    disp(TIR_GLDV_lh4x4_min);
    disp(TIR_GLDV_sd16x16);
    disp(TIR_GLDV_asm4x4_sd);
    disp(TIR_GLDV_cs4x4_max);
    disp(VIS_GLDV_con4x4_max);
    disp(TIR_GLDV_cs4x4_mean);
    disp(VIS_SADH_lh4x4_min);
    disp(VIS_RF_sre);
    disp(VIS_GLDV_ent16x16);
    disp(TIR_GLDV_asm16x16);
    disp(VIS_GLDV_ent4x4_min);
    disp(VIS_GLDV_cs4x4_max);
    disp(VIS_SF_maxpl);
    disp(VIS_SF_minpl);
    disp(TIR_SADH_cor4x4_min);
    middlequery=strcat(TIR_GLDV_lh4x4_min_val,c);
    finalquery='insert into features(TIR_GLDV_lh4x4_min_val) values(1)';
    disp(finalquery);
    curs=exec(conn,finalquery);
end