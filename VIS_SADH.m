function [mean16x16,mean4x4_max,mean4x4_min,mean4x4_mean,mean4x4_sd,sd16x16,sd4x4_max,sd4x4_min,sd4x4_mean,sd4x4_sd,asm16x16,asm4x4_max,asm4x4_min,asm4x4_mean,asm4x4_sd,con16x16,con4x4_max,con4x4_min,con4x4_mean,con4x4_sd,cor16x16,cor4x4_max,cor4x4_min,cor4x4_mean,cor4x4_sd,ent16x16,ent4x4_max,ent4x4_min,ent4x4_mean,ent4x4_sd,lh16x16,lh4x4_max,lh4x4_min,lh4x4_mean,lh4x4_sd,cs16x16,cs4x4_max,cs4x4_min,cs4x4_mean,cs4x4_sd,cp16x16,cp4x4_max,cp4x4_min,cp4x4_mean,cp4x4_sd]=VIS_SADH(vis_mat)
%caculate the 1 pixel spaced additionn matrix
[imager,imagec]=size(vis_mat);
addition_mat=ones(imager,imagec);
for addition_r=1:imager
    for addition_c=1:imagec
        if (addition_c+2<=imagec)
            addition_mat(addition_r,addition_c)=vis_mat(addition_r,addition_c)+vis_mat(addition_r,addition_c+2);
        else
            addition_mat(addition_r,addition_c)=(addition_mat(addition_r,addition_c-1)+addition_mat(addition_r,addition_c-2)+addition_mat(addition_r,addition_c-3))/3;
        end
    end
end
difference_mat=ones(imager,imagec);
for difference_r=1:imager
    for difference_c=1:imagec
        if (difference_c+2<=imagec)
            difference_mat(difference_r,difference_c)=vis_mat(difference_r,difference_c)-vis_mat(difference_r,difference_c+2);
        else
            difference_mat(difference_r,difference_c)=(difference_mat(difference_r,difference_c-1)+difference_mat(difference_r,difference_c-2)+difference_mat(difference_r,difference_c-3))/3;
        end
    end
end
%frequency matrix of each pixel value in addition matrix
ref_mat_add = unique(addition_mat);
freq_mat_add = [ref_mat_add,histc(addition_mat(:),ref_mat_add)];
tot_freq = imager*imagec;
[freq_matr_add,freq_matc_add]=size(freq_mat_add);
%frequency matrix of each pixel value in difference matrix
ref_mat_diff = unique(difference_mat);
freq_mat_diff = [ref_mat_diff,histc(difference_mat(:),ref_mat_diff)];
[freq_matr_diff,freq_matc_diff]=size(freq_mat_diff);
%calculate mean for 16x16 pixel matrix
mean16x16=0;
for m=1:freq_matr_add
    mean16x16=mean16x16+freq_mat_add(freq_matr_add,1)*(freq_mat_add(freq_matr_add,2)/tot_freq);
end

%calculate standard deviation for 16x16 pixels
sd16x16=0;
sd_16x16_add_part=0;
sd_16x16_diff_part=0;
for m=1:freq_matr_add
    sd_16x16_add_part=sd_16x16_add_part+((freq_mat_add(m,1)-mean16x16)^2)*(freq_mat_add(m,2)/tot_freq);
end
for m=1:freq_matr_diff
    sd_16x16_diff_part=sd_16x16_diff_part+((freq_mat_diff(m,1))^2)*(freq_mat_diff(m,2)/tot_freq);
end
sd16x16=sqrt((0.5)*sd_16x16_add_part+sd_16x16_diff_part);

%calculate angular second moment for 16x16 pixel matrix
asm16x16=0;
asm_16x16_add_part=0;
asm_16x16_diff_part=0;
for m=1:freq_matr_add
    asm_16x16_add_part=asm_16x16_add_part+(freq_mat_add(m,2)/tot_freq)^2;
end
for m=1:freq_matr_diff
    asm_16x16_diff_part=asm_16x16_diff_part+(freq_mat_diff(m,2)/tot_freq)^2;
end
asm16x16=asm_16x16_add_part*asm_16x16_diff_part;

%calculate contrast for 16x16 pixel matrix
con16x16=0;
for m=1:freq_matr_diff
    con16x16=con16x16+(freq_mat_diff(m,1)^2)*(freq_mat_diff(m,2)/tot_freq);
end

%calculate correlation for 16x16 pixel matrix
cor16x16=0;
cor_16x16_add_part=0;
cor_16x16_diff_part=0;
for m=1:freq_matr_add
    cor_16x16_add_part=cor_16x16_add_part+((freq_mat_add(m,1)-mean16x16)^2)*(freq_mat_add(m,2)/tot_freq);
end
for m=1:freq_matr_diff
    cor_16x16_diff_part=cor_16x16_diff_part+(freq_mat_diff(m,1)^2)*(freq_mat_diff(m,2)/tot_freq);
end
cor16x16=0.5*((cor_16x16_add_part-cor_16x16_diff_part)/((sd16x16)^2));

%calculate entropy for 16x16 pixel matrix
ent16x16=0;
ent_16x16_add_part=0;
ent_16x16_diff_part=0;
for m=1:freq_matr_add
    ent_16x16_add_part=ent_16x16_add_part+(freq_mat_add(m,2)/tot_freq)*(log10(freq_mat_add(m,2)/tot_freq));
end
for m=1:freq_matr_diff
    ent_16x16_diff_part=ent_16x16_add_part+(freq_mat_diff(m,2)/tot_freq)*(log10(freq_mat_diff(m,2)/tot_freq));
end
ent16x16=(-1)*(ent_16x16_add_part)*(ent_16x16_diff_part);

%calculate local homogeneity for 16x16 pixel matrix
lh16x16=0;
for m=1:freq_matr_diff
    lh16x16=lh16x16+(freq_mat_diff(m,2)/tot_freq)/(1+(freq_mat_diff(m,1)^2));
end

%calculate cluster shade for 16x16 pixel matrix
cs16x16=0;
for m=1:freq_matr_add
    cs16x16=cs16x16+(((freq_mat_add(m,1)-mean16x16)^3)*(freq_mat_add(m,2)/tot_freq));
end
cs16x16=(cs16x16)/((sd16x16)^3);

%calculate cluster prominence for 16x16 pixel matrix
cp16x16=0;
for m=1:freq_matr_add
    cp16x16=cp16x16+(((freq_mat_add(m,1)-mean16x16)^4)*(freq_mat_add(m,2)/tot_freq));
end
cp16x16=(cp16x16)/((sd16x16)^4)-3;
%calculate maximum, minimum, mean, standard deviation of the mean of the 16 4x4 pixel regions within the 16x16
fourxfourmat = mat2cell(vis_mat, [4,4,4,4], [4,4,4,4]);
meanvector=zeros(16,1);
sdvector=zeros(16,1);
asmvector=zeros(16,1);
convector=zeros(16,1);
corvector=zeros(16,1);
entvector=zeros(16,1);
lhvector=zeros(16,1);
csvector=zeros(16,1);
cpvector=zeros(16,1);
for a=1:4
    for b=1:4
        addition_mat_4x4=ones(4,4);
        for add_r=1:4
            for add_c=1:4
                if (add_c+2<=4)
                    addition_mat_4x4(add_r,add_c)=fourxfourmat{a,b}(add_r,add_c)+fourxfourmat{a,b}(add_r,add_c+2);
                else
                    addition_mat_4x4(add_r,add_c)=(addition_mat_4x4(add_r,add_c-1)+addition_mat_4x4(add_r,add_c-2))/2;
                end
            end
        end

        %frequency matrix of each pixel value of 4x4 matrix
        ref_mat_4x4_add = unique(addition_mat_4x4);
        freq_mat_4x4_add = [ref_mat_4x4_add,histc(addition_mat_4x4(:),ref_mat_4x4_add)];
        tot_freq_4x4=16;
        [freq_matr_4x4_add,freq_matc_4x4_add]=size(freq_mat_4x4_add);
        difference_mat_4x4=ones(4,4);
        for diff_r=1:4
            for diff_c=1:4
                if (diff_c+2<=4)
                    difference_mat_4x4(diff_r,diff_c)=fourxfourmat{a,b}(diff_r,diff_c)-fourxfourmat{a,b}(diff_r,diff_c+2);
                else
                    difference_mat_4x4(diff_r,diff_c)=(difference_mat_4x4(diff_r,diff_c-1)+difference_mat_4x4(diff_r,diff_c-2))/2;
                end
            end
        end

        %frequency matrix of each pixel value of 4x4 matrix
        ref_mat_4x4_diff = unique(addition_mat_4x4);
        freq_mat_4x4_diff = [ref_mat_4x4_diff,histc(addition_mat_4x4(:),ref_mat_4x4_diff)];
        [freq_matr_4x4_diff,freq_matc_4x4_diff]=size(freq_mat_4x4_diff);
        %calculate mean for each 4x4 matrix
        for m=1:freq_matr_4x4_add
                meanvector((a-1)*4+b,1)=meanvector((a-1)*4+b,1)+(freq_mat_4x4_add(m,1)*(freq_mat_4x4_add(m,2)/tot_freq_4x4));
        end
        
        %calculate maximum, minimum, mean, standard deviation of the standard deviation of the 16 4x4 pixel regions within the 16x16
        sd_4x4_add_part=0;
        sd_4x4_diff_part=0;
        for m=1:freq_matr_4x4_add
                sd_4x4_add_part=sd_4x4_add_part+((freq_mat_4x4_add(m,1)-mean((a-1)*4+b,1))^2)*(freq_mat_4x4_add(m,2)/tot_freq_4x4);
        end
        for m=1:freq_matr_4x4_diff
                sd_4x4_diff_part=sd_4x4_diff_part+((freq_mat_4x4_diff(m,1))^2)*(freq_mat_4x4_add(m,2)/tot_freq_4x4);
        end
        sdvector((a-1)*4+b,1)=sqrt((0.5)*sd_4x4_add_part+sd_4x4_diff_part);
        
        %calculate maximum, minimum, mean, standard deviation of the angular second moment of the 16 4x4 pixel regions within the 16x16
        asm_4x4_add_part=0;
        asm_4x4_diff_part=0;
        for m=1:freq_matr_4x4_add
                asm_4x4_add_part=asm_4x4_add_part+(freq_mat_4x4_add(m,2)/tot_freq_4x4)^2;
        end
        for m=1:freq_matr_4x4_diff
                asm_4x4_diff_part=asm_4x4_diff_part+(freq_mat_4x4_diff(m,2)/tot_freq_4x4)^2;
        end
        asmvector((a-1)*4+b,1)=asm_4x4_add_part*asm_4x4_diff_part;
       
        %calculate maximum, minimum, mean, standard deviation of the contrast of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4_diff
        convector((a-1)*4+b,1)=convector((a-1)*4+b,1)+(freq_mat_4x4_diff(m,1)^2)*(freq_mat_4x4_diff(m,2)/tot_freq_4x4);
        end
        
        %calculate maximum, minimum, mean, standard deviation of the correlation of the 16 4x4 pixel regions within the 16x16
        cor_4x4_add_part=0;
        cor_4x4_diff_part=0;
        for m=1:freq_matr_4x4_add
                cor_4x4_add_part=cor_4x4_add_part+((freq_mat_4x4_add(m,1)-mean((a-1)*4+b,1))^2)*(freq_mat_4x4_add(m,2)/tot_freq_4x4);
        end
        for m=1:freq_matr_4x4_diff
                cor_4x4_diff_part=cor_4x4_diff_part+((freq_mat_4x4_diff(m,1))^2)*(freq_mat_4x4_add(m,2)/tot_freq_4x4);
        end
        corvector((a-1)*4+b,1)=sqrt((0.5)*(cor_4x4_add_part-cor_4x4_diff_part)/(sdvector((a-1)*4+b,1)^2));
        
        %calculate maximum, minimum, mean, standard deviation of the entropy of the 16 4x4 pixel regions within the 16x16
        ent_4x4_add_part=0;
        ent_4x4_diff_part=0;
        for m=1:freq_matr_4x4_add
                ent_4x4_add_part=ent_4x4_add_part+(freq_mat_4x4_add(m,2)/tot_freq_4x4)*(log10(freq_mat_4x4_add(m,2)/tot_freq_4x4));
        end
        for m=1:freq_matr_4x4_diff
        ent_4x4_diff_part=ent_4x4_add_part+(freq_mat_4x4_diff(m,2)/tot_freq_4x4)*(log10(freq_mat_4x4_diff(m,2)/tot_freq_4x4));
        end
        entvector((a-1)*4+b,1)=(-1)*(ent_4x4_add_part)*(ent_4x4_diff_part);
        
        %calculate maximum, minimum, mean, standard deviation of the correlation of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4_diff
        lhvector((a-1)*4+b,1)=lhvector((a-1)*4+b,1)+(freq_mat_4x4_diff(m,2)/tot_freq_4x4)/(1+(freq_mat_4x4_diff(m,1)^2));
        end
        
        %calculate maximum, minimum, mean, standard deviation of the cluster shade of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4_add
        csvector((a-1)*4+b,1)=csvector((a-1)*4+b,1)+(((freq_mat_4x4_add(m,1)-meanvector((a-1)*4+b,1))^3)*(freq_mat_4x4_add(m,2)/tot_freq_4x4));
        end
        csvector((a-1)*4+b,1)=(csvector((a-1)*4+b,1))/((sdvector((a-1)*4+b,1))^3);
        
         %calculate maximum, minimum, mean, standard deviation of the cluster prominence of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4_add
        cpvector((a-1)*4+b,1)=cpvector((a-1)*4+b,1)+(((freq_mat_4x4_add(m,1)-meanvector((a-1)*4+b,1))^4)*(freq_mat_4x4_add(m,2)/tot_freq_4x4));
        end
        cpvector((a-1)*4+b,1)=(cpvector((a-1)*4+b,1))/((sdvector((a-1)*4+b,1))^4)-3;
    end
end        
mean4x4_max = max(meanvector(:));
mean4x4_min = min(meanvector(:));
mean4x4_mean = mean(double(meanvector(:)));
mean4x4_sd = std(double(meanvector(:)));
sd4x4_max = max(sdvector(:));
sd4x4_min = min(sdvector(:));
sd4x4_mean = mean(double(sdvector(:)));
sd4x4_sd = std(double(sdvector(:)));
asm4x4_max = max(asmvector(:));
asm4x4_min = min(asmvector(:));
asm4x4_mean = mean(double(asmvector(:)));
asm4x4_sd = std(double(asmvector(:)));
con4x4_max = max(convector(:));
con4x4_min = min(convector(:));
con4x4_mean = mean(double(convector(:)));
con4x4_sd = std(double(convector(:)));
cor4x4_max = max(corvector(:));
cor4x4_min = min(corvector(:));
cor4x4_mean = mean(double(corvector(:)));
cor4x4_sd = std(double(corvector(:)));
ent4x4_max = max(entvector(:));
ent4x4_min = min(entvector(:));
ent4x4_mean = mean(double(entvector(:)));
ent4x4_sd = std(double(entvector(:)));
lh4x4_max = max(lhvector(:));
lh4x4_min = min(lhvector(:));
lh4x4_mean = mean(double(lhvector(:)));
lh4x4_sd = std(double(lhvector(:)));
cs4x4_max = max(csvector(:));
cs4x4_min = min(csvector(:));
cs4x4_mean = mean(double(csvector(:)));
cs4x4_sd = std(double(csvector(:)));
cp4x4_max = max(cpvector(:));
cp4x4_min = min(cpvector(:));
cp4x4_mean = mean(double(cpvector(:)));
cp4x4_sd = std(double(cpvector(:)));
end