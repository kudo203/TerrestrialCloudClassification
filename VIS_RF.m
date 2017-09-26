function [sre,lre,gln,rln,rp]=VIS_RF(vis_mat)
[GLRLM,SI]= grayrlmatrix(vis_mat,'OFFSET',[1;2],'NumLevels',max(vis_mat(:)),'G',[]);
stats = grayrlprops(GLRLM,{'SRE','LRE','GLN','RLN','RP','LGRE','HGRE','SRLGE','SRHGE','LRLGE','LRHGE'});
sre=stats(1,1);
lre=stats(1,2);
gln=stats(1,3);
rln=stats(1,4);
rp=stats(1,5);
end