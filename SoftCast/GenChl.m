clear
clc
% SNR=[5];
tslotsnum=5;
gopcnt=16;
blk_num_i=4;
blk_num_j=3;
block_totnum=blk_num_i*blk_num_j*gopcnt;
% 
H = raylrnd(1,block_totnum,tslotsnum);
chlname=strcat('dataChl_blk',int2str(block_totnum));
save(chlname,'H');







