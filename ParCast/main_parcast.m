function main_parCast()
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

for j=1:length(SNR_table)

           for gcac=1:length(tslotTot)
             inpara.tslots_num=tslotTot(gcac);
            % SoftCast_average
            [PSNR_ParCast(j,gcac),uv,g_parcast(:,1:gcac)]=...
               mtime_latgopparfml(filename,SNR_table(j),fln,inpara.gopcnt,...
                inpara.skip,inpara.tslots_num,inpara.blk_num_i,inpara.blk_num_j,inpara.gopbg,inpara.goptmt);

            
            % SoftCast_joint
            [PSNR_ParCast_(j,gcac),uv_,g_parcast_(:,1:gcac)]=...
               mtime_latgopparfml_joint(filename,SNR_table(j),fln,inpara.gopcnt,...
                inpara.skip,inpara.tslots_num,inpara.blk_num_i,inpara.blk_num_j,inpara.gopbg,inpara.goptmt);
            
           end
                   

end

save('PSNR_ParCast','PSNR_ParCast');
save('PSNR_ParCast_','PSNR_ParCast_');

load('PSNR_ParCast');
load('PSNR_ParCast_');

figure;
plot(1:length(tslotTot),PSNR_ParCast,'rs--','LineWidth',3,'MarkerSize',12,'DisplayName','ParCast average');hold on;
plot(1:length(tslotTot),PSNR_ParCast_,'bx:','LineWidth',3,'MarkerSize',12,'DisplayName','ParCast joint');
ylabel('PSNR(dB)','FontSize',14);xlabel('Number of time slots','FontSize',14);
ylim([29 46]);
     lgd=legend('Location','se');
     lgd.FontSize=14;
     ax=gca;ax.FontSize=14;  
     ax.XAxis.TickValues=[1:1:5];                            
savefig('parcast.fig');
saveas(gcf,'parcast','epsc');


