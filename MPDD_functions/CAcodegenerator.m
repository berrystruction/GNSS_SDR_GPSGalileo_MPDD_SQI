function [time_m,codeEarly_m,codePrompt_m,codeLate_m,samplecountere,samplecounter,samplecounterl]...
    =CAcodegenerator(prn,Tc,fs,spacing,Nleft)

% spacing proportional to Tc
% work on 1ms
% PAY ATTENTION TO 0 VALUE IN COSINE LOCAL SIGNAL

tauEst=0;
Ts=1/fs;
%ns=Tc/Ts;
%fif=fs/ns;

pseudoseq=cacode(prn);
Nchip=size(pseudoseq,2); % usually 1023
fine=Nchip*Tc;
n=[0:1:floor(fine/Ts)];
tt=n*Ts;
nlength=size(n,2);
Ntot=2*Nleft+1;

% chipindex=1;
% codice=zeros(nlength,1);
% samplecounterX=zeros(Nchip,1);

% for index= 1: nlength
%     if tte(index)<=chipindex*Tc
%         codice(index)=pseudoseq(chipindex);
%         samplecounterX(chipindex)=samplecounterX(chipindex)+1;
%     else
%         chipindex=chipindex+1;
%         if chipindex>Nchip
%             chipindex=rem(chipindex,Nchip);
%         end
%         samplecounterX(chipindex)=samplecounterX(chipindex)+1;
%         codice(index)=pseudoseq(chipindex);
%     end
% end

% matrice codici shiftati
codeEarly_m=zeros(nlength,Ntot);
codeLate_m=codeEarly_m;
codePrompt_m=codeLate_m;

samplecountere=zeros(Nchip,Ntot);
samplecounter=samplecountere;
samplecounterl=samplecountere;

tt1=tt+tauEst;

time_m=zeros(3,nlength,Ntot);

for ii=-Nleft:Nleft
    
    delay=(ii)*Ts;% -Nleft-1
    ok=tt1+delay;
    
    tte=ok+spacing;
    ttp=ok;
    ttl=ok-spacing;
    
        
    time_m(1,:,ii+Nleft+1)=tte;
    time_m(2,:,ii+Nleft+1)=ttp;
    time_m(3,:,ii+Nleft+1)=ttl;


    chipindexe=1;
    chipindexp=chipindexe;
    chipindexl=chipindexp;
    chipindexl2=1;
    
%     durationE=tte(1);
%     durationP=ttp(1);
%     durationL=ttl(1);

    for index= 1: nlength     
        
        % Early
        indexCode=tte(index)/Tc;
        if indexCode<0
           indexCode=floor(indexCode);
           indexCode=rem(Nchip+indexCode+1,Nchip);
           if indexCode==0, indexCode=Nchip; end
           chipindexe=indexCode;
           codeEarly_m(index,ii+Nleft+1)=pseudoseq(chipindexe);
           samplecountere(chipindexe,ii+Nleft+1)=samplecountere(chipindexe,ii+Nleft+1)+1;

        else
            indexCode=floor(indexCode);
            indexCode=rem(Nchip+indexCode+1,Nchip);
            if indexCode==0, indexCode=Nchip; end
            chipindexe=indexCode;
            codeEarly_m(index,ii+Nleft+1)=pseudoseq(chipindexe);
            samplecountere(chipindexe,ii+Nleft+1)=samplecountere(chipindexe,ii+Nleft+1)+1;        
        end
        
        % Prompt
        indexCode=ttp(index)/Tc;
        if indexCode<0
            indexCode=floor(indexCode);
            indexCode=rem(Nchip+indexCode+1,Nchip);
            if indexCode==0, indexCode=Nchip; end
            chipindexp=indexCode;
            codePrompt_m(index,ii+Nleft+1)=pseudoseq(chipindexp);
            samplecounter(chipindexp,ii+Nleft+1)=samplecounter(chipindexp,ii+Nleft+1)+1;
            
        else
            indexCode=floor(indexCode);
            indexCode=rem(Nchip+indexCode+1,Nchip);
            if indexCode==0, indexCode=Nchip; end
            chipindexp=indexCode;
            codePrompt_m(index,ii+Nleft+1)=pseudoseq(chipindexp);
            samplecounter(chipindexp,ii+Nleft+1)=samplecounter(chipindexp,ii+Nleft+1)+1;
        end
        
        
        % Late
        indexCode=ttl(index)/Tc;
        if indexCode<0
            indexCode=floor(indexCode);
            indexCode=rem(Nchip+indexCode+1,Nchip);
            if indexCode==0, indexCode=Nchip; end
            chipindexl=indexCode;
            codeLate_m(index,ii+Nleft+1)=pseudoseq(chipindexl);
            samplecounterl(chipindexl,ii+Nleft+1)=samplecounterl(chipindexl,ii+Nleft+1)+1;
            
        else
            indexCode=floor(indexCode);
            indexCode=rem(Nchip+indexCode+1,Nchip);
            if indexCode==0, indexCode=Nchip; end
            chipindexl=indexCode;
            codeLate_m(index,ii+Nleft+1)=pseudoseq(chipindexl);
            samplecounterl(chipindexl,ii+Nleft+1)=samplecounterl(chipindexl,ii+Nleft+1)+1;
        end
        
    end   
end


% pseudoseq=[pseudoseq pseudoseq]; % elongation of the sequence
% 
% %------------------ simulazione di sinusoidi squadrate (P, E, L)
% phi=0*pi/180;% initial phase
% %sig_prompt=sign(cos(2*pi*Fcode*(tt+delay)+phi));
% delay=1*Ts;   
% sig_early=sign(cos(2*pi*Fcode*(tt-spacing+delay)+phi));
% %sig_late=sign(cos(2*pi*Fcode*(tt+spacing+delay)+phi));
% 
% init=0;
% flag=0;
% while flag==0
%     init=init+1;
%     if sig_early(init)<sig_early(init+1)
%         flag=1;
%     end
% end
% 
% %tt(init)
% n=[init-1:1:floor(fine/Ts)+init-1-1];
% tt=n*Ts;
% %tt(1)
% % SEGNALI CREATI GIUSTI
% sig_prompt=sign(cos(2*pi*Fcode*(tt+delay)+phi));
% sig_early=sign(cos(2*pi*Fcode*(tt-spacing+delay)+phi));
% sig_late=sign(cos(2*pi*Fcode*(tt+spacing+delay)+phi));
% 
% % sig_prompt(1)=1;
% % sig_early(1)=sig_prompt(1);
% % sig_late(1)=sig_early(1);
% 
% % matrice codici shiftati
% codeEarly_m=zeros(size(n,2),Ntot);
% codeLate_m=codeEarly_m;
% codePrompt_m=codeLate_m;
% 
% codeEarly_m(:,1)=sig_early;
% codePrompt_m(:,1)=sig_prompt;
% codeLate_m(:,1)=sig_late;
% 
% sig_early_m(:,1)=sig_early;
% sig_prompt_m(:,1)=sig_prompt;
% sig_late_m(:,1)=sig_late;
% 
% 
% for ii=2:Ntot
%     
%     delay=ii*Ts;  
%     
%     sig_prompt=sign(cos(2*pi*Fcode*(tt+delay)+phi));
%     sig_early=sign(cos(2*pi*Fcode*(tt-spacing+delay)+phi));
%     sig_late=sign(cos(2*pi*Fcode*(tt+spacing+delay)+phi));
% 
%     codeEarly_m(:,ii)=sig_early;
%     codePrompt_m(:,ii)=sig_prompt;
%     codeLate_m(:,ii)=sig_late;
%     
%     sig_early_m(:,ii)=sig_early;
%     sig_prompt_m(:,ii)=sig_prompt;
%     sig_late_m(:,ii)=sig_late;
% 
% end
% 
% 
% 
% 
% for ii=1:Ntot
%     
%     % figure
%     % plot(tt,sig_early,'r^-',tt,sig_prompt,'bo-',tt,sig_late,'g*-')
%     % grid on
%     % legend('sig_early','sig_prompt','sig_late');
%     
%     seqindex=1;
%     seqindexe=seqindex;
%     seqindexl=seqindex;
%     
%     codeprompt=zeros(size(sig_prompt));
%     codeearly=codeprompt;
%     codelate=codeprompt;
%     
%     samplecounter=zeros(Nchip,1);
%     samplecountere=samplecounter;
%     samplecounterl=samplecounter;
%     
%     
%     for indexchip=1:length(n)-1
%         % Prompt
%         if seqindex>Nchip
%             seqindex=rem(seqindex,Nchip);
%         end
%         if codePrompt_m(indexchip,ii)>=codePrompt_m(indexchip+1,ii)
%             codeprompt(indexchip)=pseudoseq(seqindex);
%             %samplecounter(seqindex)=samplecounter(seqindex)+1;
%         else
%             codeprompt(indexchip)=pseudoseq(seqindex);
%             %samplecounter(seqindex)=samplecounter(seqindex)+1;
%             seqindex=seqindex+1;
%         end
%         
%         % Early
%         if seqindexe>Nchip
%             seqindexe=rem(seqindexe,Nchip);
%         end
%         if codeEarly_m(indexchip,ii)>=codeEarly_m(indexchip+1,ii)
%             codeearly(indexchip)=pseudoseq(seqindexe);
%             samplecountere(seqindexe)=samplecountere(seqindexe)+1;
%         else
%             codeearly(indexchip)=pseudoseq(seqindexe);
%             %samplecountere(seqindexe)=samplecountere(seqindexe)+1;
%             seqindexe=seqindexe+1;
%         end
%         
%         % Late
%         if seqindexl>Nchip
%             seqindexl=rem(seqindexl,Nchip);
%         end
%         
%         if codeLate_m(indexchip,ii)>=codeLate_m(indexchip+1,ii)
%             codelate(indexchip)=pseudoseq(seqindexl);
%             %samplecounterl(seqindexl)=samplecounterl(seqindexl)+1;
%         else
%             codelate(indexchip)=pseudoseq(seqindexl);
%             %samplecounterl(seqindexl)=samplecounterl(seqindexl)+1;
%             seqindexl=seqindexl+1;
%         end
%         
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     %indexchip=indexchip+1;
%     % Prompt
%     if codePrompt_m(indexchip,ii)<=codePrompt_m(indexchip+1,ii)
%         codeprompt(indexchip)=pseudoseq(seqindex);
%         %samplecounter(seqindex)=samplecounter(seqindex)+1;
%     else
%         codeprompt(indexchip)=pseudoseq(1);
%        % samplecounter(1)=samplecounter(1)+1;
%     end
%     
%     % Early
%     if codeEarly_m(indexchip,ii)<=codeEarly_m(indexchip+1,ii)
%         codeearly(indexchip)=pseudoseq(seqindexe);
%         %samplecountere(seqindexe)=samplecountere(seqindexe)+1;
%     else
%         codeearly(indexchip)=pseudoseq(1);
%         %samplecountere(1)=samplecountere(1)+1;
%     end
%     
%     % Late
%     if codeLate_m(indexchip,ii)<=codeLate_m(indexchip+1,ii)
%         codelate(indexchip)=pseudoseq(seqindexl);
%         %samplecounterl(seqindexl)=samplecounterl(seqindexl)+1;
%     else
%         codelate(indexchip)=pseudoseq(1);
%         %samplecounterl(1)=samplecounterl(1)+1;
%     end
%   
%     
%     codeEarly_m(:,ii)=codeearly;
%     codePrompt_m(:,ii)=codeprompt;
%     codeLate_m(:,ii)=codelate;
% 
% end
% 
% figure
% dd=10;
% subplot(311)
% plot(n,codeEarly_m(:,cc),'b.-',n,codice,'ro-')
% title(['sequenza chip originali ' num2str(pseudoseq(1:dd))])
% legend('codeEarly','codice')
% grid on
% 
% subplot(312)
% plot(n,codePrompt_m(:,cc),'b.-',n,codice,'ro-')
% title(['sequenza chip originali ' num2str(pseudoseq(1:dd))])
% legend('codePrompt','codice')
% grid on
% 
% subplot(313)
% plot(n,codeLate_m(:,cc),'b.-',n,codice,'ro-')
% title(['sequenza chip originali ' num2str(pseudoseq(1:dd))])
% legend('codeLate','codice')
% grid on

% figure
% plot([1:Nchip],samplecounterX)
% title(['# of samples per chip. Sequenza chip originali ' num2str(pseudoseq(1:dd))])
% grid on
end