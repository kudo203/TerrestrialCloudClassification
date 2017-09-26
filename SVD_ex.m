function a=SVD_ex()
a=1;
VIS_U=zeros(8,8);
VIS_S=zeros(8,8);
VIS_V=zeros(8,8);

TIR_U=zeros(8,8);
TIR_S=zeros(8,8);
TIR_V=zeros(8,8);

S=zeros(8,8);
data_VIS=h5read('K1VHR_01AUG2010_1130_L02_ASI.h5','/VHRR/Image Data/VHRR_VIS');
data_TIR=h5read('K1VHR_01AUG2010_1130_L02_ASI.h5','/VHRR/Image Data/VHRR_TIR');

matrix_VIS = mat2cell(data_VIS, [512,297], [512,296]);
matrix1_VIS= mat2cell(matrix_VIS{1,1},[8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8],[8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]);

matrix_TIR = mat2cell(data_TIR, [512,297], [512,296]);
matrix1_TIR= mat2cell(matrix_TIR{1,1},[8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8],[8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]);

for m=1:8
    for n=1:8
        [VIS_U,VIS_S,VIS_V]=svd(double(matrix1_VIS{m,n}));
        [TIR_U,TIR_S,TIR_V]=svd(double(matrix1_TIR{m,n}));
        for a=1:8
        VIS_S(a,a)=log10(VIS_S(a,a));
        TIR_S(a,a)=log10(TIR_S(a,a));
        S(a,a)=(VIS_S(a,a)+TIR_S(a,a))/2;
        end
        disp(S);
        fprintf('\n');
    end
end
end