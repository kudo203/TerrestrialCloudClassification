function [MAX,MIN,MEAN,MEDIAN,MODE,RANGE,SD]=VIS_SF(vis_mat)
MAX = max(vis_mat(:));
MIN = min(vis_mat(:));
MEAN=mean(double(vis_mat(:)));
MEDIAN=median(double(vis_mat(:)));
MODE=mode(double(vis_mat(:)));
RANGE=MAX-MIN;
SD=std(double(vis_mat(:)));
end 