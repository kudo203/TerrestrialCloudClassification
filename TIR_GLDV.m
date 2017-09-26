function [mean16x16,mean4x4_max,mean4x4_min,mean4x4_mean,mean4x4_sd,sd16x16,sd4x4_max,sd4x4_min,sd4x4_mean,sd4x4_sd,asm16x16,asm4x4_max,asm4x4_min,asm4x4_mean,asm4x4_sd,ent16x16,ent4x4_max,ent4x4_min,ent4x4_mean,ent4x4_sd,lh16x16,lh4x4_max,lh4x4_min,lh4x4_mean,lh4x4_sd,con16x16,con4x4_max,con4x4_min,con4x4_mean,con4x4_sd,cs16x16,cs4x4_max,cs4x4_min,cs4x4_mean,cs4x4_sd,cp16x16,cp4x4_max,cp4x4_min,cp4x4_mean,cp4x4_sd]=TIR_GLDV(ir_mat)
%caculate the 1 pixel spaced difference matrix
[imager,imagec]=size(ir_mat);
difference_mat=ones(imager,imagec);
for difference_r=1:imager
    for difference_c=1:imagec
        if (difference_c+2<=imagec)
            difference_mat(difference_r,difference_c)=abs(ir_mat(difference_r,difference_c)-ir_mat(difference_r,difference_c+2));
        else
            difference_mat(difference_r,difference_c)=(difference_mat(difference_r,difference_c-1)+difference_mat(difference_r,difference_c-2)+difference_mat(difference_r,difference_c-3))/3;
        end
    end
end
%frequency matrix of each pixel value
ref_mat = unique(difference_mat);
freq_mat = [ref_mat,histc(difference_mat(:),ref_mat)];
tot_freq=imager*imagec;
%calculate mean for 16x16 pixels
mean16x16=0;
[freq_matr,freq_matc]=size(freq_mat);
for m=1:freq_matr
        mean16x16=mean16x16+(freq_mat(m,1)*(freq_mat(m,2)/tot_freq));
end

%calculate standard deviation for 16x16 pixels
sd16x16=0;
for m=1:freq_matr
        sd16x16=sd16x16+(((freq_mat(m,1)-mean16x16)^2)*(freq_mat(m,2)/tot_freq));
end
sd16x16=sqrt(sd16x16);

%calculate angular second moment for 16x16 pixels
asm16x16=0;
for m=1:freq_matr
        asm16x16=asm16x16+((freq_mat(m,2)/tot_freq)^2);
end

%calculate entropy for 16x16 pixels
ent16x16=0;
for m=1:freq_matr
        ent16x16=ent16x16+((freq_mat(m,2))/tot_freq)*(log10((freq_mat(m,2))/tot_freq));
end
ent16x16=(-1)*ent16x16;

%calculate local homogeneity for 16x16 pixels
lh16x16=0;
for m=1:freq_matr
        lh16x16=lh16x16+(freq_mat(m,2)/tot_freq)/(1+(freq_mat(m,1))^2);
end

%calculate contrast for 16x16 pixels
con16x16=0;
for m=1:freq_matr
        con16x16=con16x16+((freq_mat(m,1)^2)*(freq_mat(m,2)/tot_freq));
end

%calculate cluster shade for 16x16 pixels
cs16x16=0;
for m=1:freq_matr
         cs16x16=cs16x16+(((freq_mat(m,1)-mean16x16)^3)*(freq_mat(m,2)/tot_freq));
end
cs16x16=cs16x16/(sd16x16^3);

%calculate cluster prominence for 16x16 pixels
cp16x16=0;
for m=1:freq_matr
         cp16x16=cp16x16+(((freq_mat(m,1)-mean16x16)^4)*(freq_mat(m,2)/tot_freq));
end
cp16x16=(cp16x16/(sd16x16^4))-3;
%calculate maximum, minimum, mean, standard deviation of the mean of the 16 4x4 pixel regions within the 16x16
fourxfourmat = mat2cell(ir_mat, [4,4,4,4], [4,4,4,4]);
meanvector=zeros(16,1);
sdvector=zeros(16,1);
asmvector=zeros(16,1);
entvector=zeros(16,1);
lhvector=zeros(16,1);
convector=zeros(16,1);
csvector=zeros(16,1);
cpvector=zeros(16,1);
for a=1:4
    for b=1:4
        difference_mat_4x4=ones(4,4);
        for diff_r=1:4
            for diff_c=1:4
                if (diff_c+2<=4)
                    difference_mat_4x4(diff_r,diff_c)=abs(fourxfourmat{a,b}(diff_r,diff_c)-fourxfourmat{a,b}(diff_r,diff_c+2));
                else
                    difference_mat_4x4(diff_r,diff_c)=(difference_mat_4x4(diff_r,diff_c-1)+difference_mat_4x4(diff_r,diff_c-2))/2;
                end
            end
        end

        %frequency matrix of each pixel value of 4x4 matrix
        ref_mat_4x4 = unique(difference_mat_4x4);
        freq_mat_4x4 = [ref_mat_4x4,histc(difference_mat_4x4(:),ref_mat_4x4)];
        tot_freq_4x4=16;
        %calculate mean for each 4x4 matrix
        [freq_matr_4x4,freq_matc_4x4]=size(freq_mat_4x4);
        for m=1:freq_matr_4x4
                meanvector((a-1)*4+b,1)=meanvector((a-1)*4+b,1)+(freq_mat_4x4(m,1)*(freq_mat_4x4(m,2)/tot_freq_4x4));
        end
        %calculate maximum, minimum, mean, standard deviation of the standard deviation of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4
                sdvector((a-1)*4+b,1)=sdvector((a-1)*4+b,1)+(((freq_mat_4x4(m,1)-meanvector((a-1)*4+b,1))^2)*(freq_mat_4x4(m,2)/tot_freq_4x4));
        end
        sdvector((a-1)*4+b,1)=sqrt(sdvector((a-1)*4+b,1));
        
        %calculate maximum, minimum, mean, standard deviation of the angular second moment of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4
                asmvector((a-1)*4+b,1)=asmvector((a-1)*4+b,1)+((freq_mat_4x4(m,2)/tot_freq_4x4)^2);
        end
        
        %calculate maximum, minimum, mean, standard deviation of the entropy of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4
                entvector((a-1)*4+b,1)=entvector((a-1)*4+b,1)+((freq_mat_4x4(m,2))/tot_freq_4x4)*(log10((freq_mat_4x4(m,2))/tot_freq_4x4));
        end
        entvector((a-1)*4+b,1)=(-1)*entvector((a-1)*4+b,1);
    
        %calculate maximum, minimum, mean, standard deviation of the entropy of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4
                lhvector((a-1)*4+b,1)=lhvector((a-1)*4+b,1)+(freq_mat_4x4(m,2)/tot_freq_4x4)/(1+(freq_mat_4x4(m,1))^2);
        end
    
        %calculate maximum, minimum, mean, standard deviation of the contrast of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4
                convector((a-1)*4+b,1)=convector((a-1)*4+b,1)+((freq_mat_4x4(m,1)^2)*(freq_mat_4x4(m,2)/tot_freq_4x4));
        end
        %calculate maximum, minimum, mean, standard deviation of the cluster shade of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4
                csvector((a-1)*4+b,1)=csvector((a-1)*4+b,1)+(((freq_mat_4x4(m,1)-meanvector((a-1)*4+b,1))^3)*(freq_mat_4x4(m,2)/tot_freq_4x4));
        end
        csvector((a-1)*4+b,1)=csvector((a-1)*4+b,1)/(sdvector((a-1)*4+b,1)^3);    
    

        %calculate maximum, minimum, mean, standard deviation of the cluster prominenece of the 16 4x4 pixel regions within the 16x16
        for m=1:freq_matr_4x4
                cpvector((a-1)*4+b,1)=cpvector((a-1)*4+b,1)+(((freq_mat_4x4(m,1)-meanvector((a-1)*4+b,1))^4)*(freq_mat_4x4(2,2)/tot_freq_4x4));
        end
        cpvector((a-1)*4+b,1)=(cpvector((a-1)*4+b,1)/(sdvector((a-1)*4+b,1)^4))-3;
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
        
ent4x4_max = max(entvector(:));
ent4x4_min = min(entvector(:));
ent4x4_mean = mean(double(entvector(:)));
ent4x4_sd = std(double(entvector(:)));

lh4x4_max = max(lhvector(:));
lh4x4_min = min(lhvector(:));
lh4x4_mean = mean(double(lhvector(:)));
lh4x4_sd = std(double(lhvector(:)));

con4x4_max = max(convector(:));
con4x4_min = min(convector(:));
con4x4_mean = mean(double(convector(:)));
con4x4_sd = std(double(convector(:)));

cs4x4_max = max(csvector(:));
cs4x4_min = min(csvector(:));
cs4x4_mean = mean(double(csvector(:)));
cs4x4_sd = std(double(csvector(:)));

cp4x4_max = max(cpvector(:));
cp4x4_min = min(cpvector(:));
cp4x4_mean = mean(double(cpvector(:)));
cp4x4_sd = std(double(cpvector(:)));

end