function main_softCast()
clear all
clc

%%   parameter setting
%filelist={'/Users/saqr/Documents/VIDEOS/foreman_cif/foreman_cif.yuv'};
fln='fm';           % stands for what that is related to AWGN
SNR_table=[5];      % SNR=[-5,0,5,10]
inpara.gopcnt=16;        % number of Framesin 1 GOP (central GOP)  % number of frames in GOP
inpara.blk_num_i=4;      % number of chunks Horizontally so Chunk width=88
inpara.blk_num_j=3;      % number of chunks Vertically Chunk height=96 not 72
inpara.skip=1;
inpara.gopbg=1;          % GOP of beginning or begin frame of GOP   
inpara.goptmt=inpara.gopbg+inpara.gopcnt-1; % GOP of terminating

n_width=352;
n_height=288;
inpara.n_width=352;
inpara.n_height=288;

blk_size_i = n_width/inpara.blk_num_i;
blk_size_j = n_height/inpara.blk_num_j;
block_size=blk_size_i*blk_size_j;
picsize = n_width*n_height;

tslotTot=[1,2,3,4,5];   % time slot number=1 MCast primary property
filename='/Users/saqr/Documents/VIDEOS/foreman_cif/foreman_cif.yuv';
fid = fopen(filename,'r');
%outname=strcat('outsnr', int2str(SNR_table) ,'ts',int2str(tslotTot),'_Mcast.yuv');

blktotnum=inpara.gopcnt*inpara.blk_num_i*inpara.blk_num_j;
blktotidx=[1:blktotnum];

% RayLeigh Fading 
ldn=strcat('dataChl_blk',int2str(blktotnum)); %ldn=strcat('rayleighb1_',int2str(block_totnum));
load(ldn);
H=H;   % The overall channel materix of all subcarriers stacked in one

for frame =1:inpara.gopcnt
    fseek(fid, (1+frame-2)*inpara.skip*picsize*1.5, -1);  
    I = fread (fid, n_width * n_height, 'uint8');
    I = reshape(I, n_width, n_height);
    I3D(:,:,frame)=I;
    imshow(uint8(I3D(:,:,frame)));
end
show_frames_TD(I3D,inpara.gopcnt,1);    %SAQR 
sigma_signal = (norm(I3D(:))^2/(picsize*inpara.gopcnt));              % std of signal


for j=1:length(SNR_table)
    sigma_noise = sqrt( sigma_signal/10^(0.1 * SNR_table(j)) );  % std of noise
    powertot=sigma_signal*blktotnum;
    noise=gnoise(sigma_noise,SNR_table(j),inpara);

           for gcac=1:length(tslotTot)
             inpara.tslots_num=tslotTot(gcac);
            % SoftCast_average
            [PSNR_SoftCast(j,gcac),noise,uv,g_softcast(:,1:gcac)]=...
                mtime_latgopfml(filename,noise,H,SNR_table(j),fln,inpara.gopcnt,...
                inpara.skip,inpara.tslots_num,inpara.blk_num_i,inpara.blk_num_j,inpara.gopbg,inpara.goptmt);
           % SoftCast_joint
            [PSNR_SoftCast_(j,gcac),noise_,uv_,g_softcast_(:,1:gcac)]=...
                mtime_latgopfml_joint(filename,noise,H,SNR_table(j),fln,inpara.gopcnt,...
                inpara.skip,inpara.tslots_num,inpara.blk_num_i,inpara.blk_num_j,inpara.gopbg,inpara.goptmt);
            
           end
                   

end

save('PSNR_SoftCast','PSNR_SoftCast');
save('PSNR_SoftCast_','PSNR_SoftCast_');

load('PSNR_SoftCast');
load('PSNR_SoftCast_');

figure;
plot(1:length(tslotTot),PSNR_SoftCast_,'g.:','LineWidth',3,'MarkerSize',38,'DisplayName','SoftCast joint');hold on;
plot(1:length(tslotTot),PSNR_SoftCast,'bo--','LineWidth',3,'MarkerSize',12,'DisplayName','SoftCast average');
ylabel('PSNR(dB)','FontSize',14);xlabel('Number of time slots','FontSize',14);
ylim([29 46]);
     lgd=legend('Location','se');
     lgd.FontSize=14;
     ax=gca;ax.FontSize=14;
     ax.XAxis.TickValues=[1:1:5];  
savefig('softcast.fig');
saveas(gcf,'softcast','epsc');


