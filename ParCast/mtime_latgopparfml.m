function [PSNRgop,uv,g]=mtime_latgopparfml(filename,SNR_table,fln,gopcnt,skip,tslots_num,blk_num_i,blk_num_j,gopbg,goptmt)
%% parameter setting
SNR_list=SNR_table;
block_num=blk_num_i*blk_num_j;     % total number of blocks in 1 frame 
block_totnum=block_num*gopcnt;     % total number of blocks in 1 GOP
txbands=block_totnum;              % Transmission Bands
txtotbands=block_totnum*tslots_num;% This Data will be sent only once

% RayLeigh Fading 
ldn=strcat('dataChl_blk',int2str(block_totnum)); %ldn=strcat('rayleighb1_',int2str(block_totnum));
load(ldn);
H=H;   % The overall channel materix of all subcarriers stacked in one

% AWGN
gnname=strcat('gnsnr',int2str(SNR_table(1)),'_blk',int2str(txbands),'_',fln);
load(gnname);
gnoise=gnoise;  % Gaussian Noise AWGN for 48 channels/subcarriers

gop_nobg=gopbg;
gop_no=gopbg;
goptmt=goptmt;

hbidx= ([1:block_totnum]');
SNR_list = SNR_table;
PSNR = zeros(length(SNR_list), gopcnt);

outname=strcat('outsnr', int2str(SNR_table) ,'ts',int2str(tslots_num),'_parcast.yuv');
%% Digital Modulation
n_width=352;w=352;        % Video Width 
n_height=288;h=288;         % Video Height
picsize = n_width*n_height;        % Number of coefficient in 1 Frame 

blk_size_i = n_width/blk_num_i;      % Chunk width=88
blk_size_j = n_height/blk_num_j;     % Chunk height=96 288/3=72
block_size=blk_size_i*blk_size_j;    % Number of coefficient in 1 chunk 

for hhcac=1:tslots_num                % Prepare channel Gain based on time slots
    hh=H(:,hhcac);       % Power allocation by treating H as parameters
    htx(:,hhcac)=sort(hh,'descend');         % Order Channel Gain 
end 

% test image = 32x32 piece of cameraman's arm
fid = fopen(filename,'r');
fid_out = fopen(outname,'w');
fseek(fid_out, 0, -1);  
q_mod_array = zeros(length(SNR_list),gopcnt);

I3D=zeros(n_width,n_height,gopcnt);
   
for iSNR = 1:length(SNR_list)
        iSNR 
    for gop=gop_nobg:16:goptmt 
        uv=uint8(zeros(n_width*n_height/2,gopcnt));
        for frame =1:gopcnt
            fseek(fid, (gop+frame-2)*skip*picsize*1.5, -1);  
            I = fread (fid, n_width * n_height, 'uint8');
            uv(:,frame) = fread (fid, picsize/2, 'uint8');
            I = reshape(I, n_width, n_height);
            I3D(:,:,frame)=I;
            imshow(uint8(I3D(:,:,frame)));
        end
        show_frames_TD(I3D,gopcnt,1);    %SAQR 
        sigma_signal = (norm(I3D(:))^2/(picsize*gopcnt));              % std of signal
        sigma_noise = sqrt( sigma_signal/10^(0.1 * SNR_list(iSNR)) );  % std of noise
        power_noise=sigma_noise^2;                                     % variance of noise

        powertot=sigma_signal*block_totnum; % Total power = std of signal x total number of chunks   
        % powertot=powertot*tslots_num;
        
        I_raw = reshape(I3D, n_width * n_height, gopcnt); % conversion from 3D to stack of 1D. (each frame 1D)
        I3D=I3D-128;                       % Level-offsetting to Get a zero centered Data
        
        %% 2D-DCT Spatial Decorrelation
        for frame =1:gopcnt
            I3D(:,:,frame)=dct(dct(I3D(:,:,frame))')';         
        end
       %% 1D-DCT Temporal Decorrelation
        for i=1:w
            for j=1:h
                I3D(i,j,:)=dct(I3D(i,j,:));
            end
        end
        show_frames_FD(I3D,gopcnt,2);
        
     %% 3D to 2D Conversion   
        lumbda = zeros(block_totnum, 1);
        I_large=[];
        for gfcac=1:gopcnt                          
          I_large=[I_large,I3D(:,:,gfcac)]; %Each frame is next to each other, all have same height
        end
        x=I_large; % Heightx(Width*GOP)
     %% Chunk Division
     % wierd method to calculate the variance of the blocks but it is still accurate, I examined it.
     % why wierd? because each coefficient location for all chunks is saved in 1 chunk, e.g., 1st coeff. for all chunks is saved in 1st chunk
     % so now each chunk is not gathering the neighboring coeff. but chunks size is 4x12 instead of 88x96. But eventually it is similar to mine
        for ii=1:blk_size_i
            for jj=1:blk_size_j
                x_blk=x(ii:blk_size_i:n_width, jj:blk_size_j:n_height*gopcnt);
                x_vec=reshape( x_blk,block_totnum,1);
                lumbda = lumbda + x_vec.^2;
            end
        end
        lumbda = lumbda/(blk_size_i*blk_size_j);
        % all of this weird calculation but it is still yeild to the same result like mine.
        [lumbda_sorted, sortblk_idx]=sort(lumbda,'descend');  % Order Chunk variance

        %% ParCast
        g=zeros(block_totnum,tslots_num); % to accumulate all bij
        for gcac=1:tslots_num     % the iteration depends how many time slots we have 
                                  % every place we see gcac it is a sign for bij 
              htxtt=htx(:,gcac);  % we only have 1 48-channel gain stock already reorder
               % Scaling factor
              gfml=gpafml(htxtt,lumbda_sorted,powertot,block_totnum,hbidx);
              g(:,gcac)=gfml;

              I_rec1 = zeros(n_width, n_height*gopcnt);    
                 for bb=1:block_totnum
                     curblk_index=sortblk_idx(bb); % index of real locations of the reorderred chunks
                     curblk_index=curblk_index-1;     
                     m_divider=floor(curblk_index/blk_num_i);
                     n_rmnd=mod(curblk_index,blk_num_i);
                     loc_i=n_rmnd*blk_size_i+1;      % used to locate the exact location of chunks,
                     loc_j=m_divider*blk_size_j+1;   % after being reordered to grap its content

                     sigma_scur=lumbda_sorted(bb); % lumbdai
                     loc_blkb=find(hbidx==bb); % decides what G,H,AWGN I selected for each specific chunk
                     if sum(loc_blkb)<0
                           I_rec1(loc_i:(loc_i+blk_size_i-1), loc_j:(loc_j+blk_size_j-1))= 0;
                           continue
                     end       
                     gcur=gfml(loc_blkb);     % Gi
                     htxcur= htxtt(loc_blkb); % Hi

                     Hgi= htxcur.*gcur; %H*G*R                
                     x_bblk=x(loc_i:(loc_i+blk_size_i-1), loc_j:(loc_j+blk_size_j-1)); % Chunk content
                     xx_blk=x_bblk(:)'; % 2D to 1D
                     m = xx_blk;        % X only 1D vector
                     %e1=sigma_noise*gnoise(loc_blkb.*gcac,:);   %SAQR
                     
                     e1=gnoise(loc_blkb.*gcac,:); %When tslots_num=2, AWGN become just the even lines which is not sufficient, w
                     y1 = Hgi*m + e1;             % H*G*R*X+AWGN    the received signal at the first receive antenna
                                                  % Yi=(Hi*Gi)*Xi+Zi MCast equation(2)
                     %% MMSE 
                     w1=sum(Hgi.^2); %|H.G|^2 The Sum is for establishing bij=sum(Hgi(loc_blkb,:) w1=sum(Hgi(loc_blkb,:).^2);
                     w2=w1*sigma_scur+sigma_noise^2; %lambda*|H.G|^2+n_std^2
                     mmse_w=sigma_scur*(Hgi)'/w2;                                 % MCast equ(3)without Yi
                     mse_mlatheory(iSNR,bb,gop_no)=sigma_noise^2/(w2/sigma_scur); % MCast equ(4)
                     xx_dec=mmse_w*y1;                                            % MCast equ(3) complete 
                     xx_dec=xx_dec';                              
                     x_dec=reshape(xx_dec,blk_size_i,blk_size_j); %  xx_blk=x_bblk(:)'; % 2D to 1D
                     mse_mlatprac(iSNR,bb,gop_no)=sum((xx_dec'-xx_blk).^2)/(blk_size_i*blk_size_j); % compare MSE for each chunk in SNR
                     I_rec1(loc_i:(loc_i+blk_size_i-1), loc_j:(loc_j+blk_size_j-1))= x_dec;         % rebuild the 2D matrix       
                 end   
                 
                 %% 2D to 3D Conversion
                    for gprc=1:gopcnt
                        I3D(:,:,gprc)=I_rec1(1:n_width,n_height*(gprc-1)+1:n_height*gprc); 
                    end  
                    
                 %% 3D to 4D Conversion
                    I3Drcd(:,:,:,gcac)=I3D; 

        end    %for every additional time slot

                                         % ParCast+Joint 
                                         % the receiver jointly decodes the received signals of all time slots.
         I3Dfl=sum(I3Drcd,4)./tslots_num;% ParCast+Average
                                         % the receiver decodes the received signal of each time slot independently and average the results
         %% 1D-iDCT Temporal Decorrelation
                    for i=1:w
                        for j=1:h
                            I3Dfl(i,j,:)=idct(I3Dfl(i,j,:));
                        end
                    end
                    
        %% 2D-iDCT Spatial Decorrelation
                    for frame =1:gopcnt
                        I3Dfl(:,:,frame)=idct(idct(I3Dfl(:,:,frame)')');
                    end

                    I3Dfl = I3Dfl + 128;
                    I3Dfl = uint8(I3Dfl);                 
                    show_frames_TD(I3Dfl,gopcnt,4);    %SAQR 
                    I3Dfl = double(I3Dfl);
                    I_rec1 = reshape( I3Dfl, n_width * n_height, gopcnt);            
                    
                    for frame=1:gopcnt
                         MSE(iSNR, gop+frame-1) =(norm(I_rec1(:,frame)-I_raw(:,frame))^2/length(I_rec1));
                         PSNR(iSNR, gop+frame-1) = 10 * log10(255^2/MSE(iSNR, gop+frame-1));
                         PSNR(iSNR, gop+frame-1);    
                         fwrite(fid_out, I_rec1(:,frame), 'uint8');
                         fwrite(fid_out, uv(:,frame) , 'uint8');   
                    end
                    fclose(fid); %SAQR 
                    fclose(fid_out);
                    
                %PSNRgop(iSNR, gop_no)=mean(PSNR);
                PSNRgop(iSNR, gop_no)=sum(PSNR(iSNR,:))/gopcnt; 
                PSNRgop(iSNR, gop_no)

    end%gop


end%SNR
%fclose(fid);     SAQR commented this
%fclose(fid_out);

% return
%rst=PSNRgop;  %SAQR commented this 
end





                                                            
