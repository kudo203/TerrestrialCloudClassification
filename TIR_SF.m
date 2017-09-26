function [MAX,MIN,MEAN,MEDIAN,MODE,RANGE,SD]=TIR_SF(ir_mat)
MAX = max(ir_mat(:));
MIN = min(ir_mat(:));
MEAN=mean(double(ir_mat(:)));
MEDIAN=median(double(ir_mat(:)));
MODE=mode(double(ir_mat(:)));
RANGE=MAX-MIN;
SD=std(double(ir_mat(:)));
end 
 
 

 