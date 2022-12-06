function gnoise=gnoise(sigma_noise,SNR,inpara)
n_width=inpara.n_width;
n_height=inpara.n_height;
gopcnt=inpara.gopcnt;
blk_num_i=inpara.blk_num_i;
blk_num_j=inpara.blk_num_j;

tslotsnum=5;
block_totnum=blk_num_i*blk_num_j*gopcnt;

blk_size_i = n_width/blk_num_i;      % Chunk width=88
blk_size_j = n_height/blk_num_j;     % Chunk height=96 288/3=72
block_size=blk_size_i*blk_size_j;    % Number of coefficient in 1 chunk 
% 

gnoise=sigma_noise*randn(block_totnum*tslotsnum,block_size);

text=strcat('gnsnr',int2str(SNR),'_blk',int2str(block_totnum),'_','fm');
save(text,'gnoise');
end 







