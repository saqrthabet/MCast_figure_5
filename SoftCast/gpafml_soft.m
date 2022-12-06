function gfml=gpafml_soft(htxtt,lumbda,powertot,block_totnum,hbidx)
% addpath('C:\Users\HCF\Desktop\sfcode_opt8blk64');
% addpath('C:\Users\HCF\Desktop\sfcode_opt8blk64\samples');

gg=zeros(block_totnum,1);              
wg2=0;
     % Sum( sqrt(lambda)/s )
     for chcpa=1:block_totnum                     
         locblkpa=hbidx(chcpa);
         if locblkpa<1
                 continue
         end 
         sigma_sgpa=lumbda(locblkpa);   % single Chunk variance
         sigma_sgpa=sqrt(sigma_sgpa);          % single sqrt Chunk variance
         hcurpa=htxtt(chcpa);                  % single Channel Gain
         wg2=wg2+sigma_sgpa;%/hcurpa;            % sqrt(lambda)/s
     end 
     wg2=wg2;
    % the rest of equation 3 in ParCast 
     for chc=1:block_totnum
           locblkchc=hbidx(chc);
           if locblkchc<1
                gg(chc)=0;
                continue
           end 
       sigma_sg=lumbda(locblkchc);
       sigma_sg=sqrt(sigma_sg);
       hcurchc=htxtt(chc);
       wg1=sigma_sg;%*hcurchc;         % sqrt(lambda)*s
       gmda=sqrt(powertot/wg1/wg2);
       gg(chc)=gmda;                
     end
     
     gfml=gg;
     g=gg;

end