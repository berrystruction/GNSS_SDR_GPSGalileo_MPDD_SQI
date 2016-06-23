% function results=countSuccessDetection(indices,mvgAvgtime,mvgAvgtime2,tracktime,skip_sec)
% Npoints=tracktime/mvgAvgtime/mvgAvgtime2;
% counterMP=zeros(1,5);
% sogliaMP1=(15-skip_sec)/(mvgAvgtime*1e-3)/mvgAvgtime2; % indice per sapere qundo inizia il MP neli dataset NFUELS
% sogliaMP2=(35-skip_sec)/(mvgAvgtime*1e-3)/mvgAvgtime2; %
% durationMP=5/(mvgAvgtime*1e-3); % duration of the MP (5s) in number of points

function results=countSuccessDetection(indices,mvgAvgtime,mvgAvgtime2,timesize,skip_sec)
Npoints=timesize;
counterMP=zeros(1,5);
sogliaMP1=(15-skip_sec)/(mvgAvgtime*1e-3)/mvgAvgtime2; % indice per sapere qundo inizia il MP neli dataset NFUELS
sogliaMP2=(35-skip_sec)/(mvgAvgtime*1e-3)/mvgAvgtime2; %
durationMP=5/(mvgAvgtime*1e-3); % duration of the MP (5s) in number of points
LOSindex=1; % theoretical index vector of the signal with only LOS

if sogliaMP2>Npoints
    sogliaMP2=sogliaMP1;
end

for index=1:sogliaMP1-1
    if indices(index)~=LOSindex
        counterMP(1)=counterMP(1)+1;
    end
end
disp('-----------------------------------------------------------------------------------------------')
disp(['Ho contato ' num2str(counterMP(1)) ' MP points su ' num2str(sogliaMP1) ' totali (only noise)'])

%%%
for index=sogliaMP1:sogliaMP1+durationMP-1
    if indices(index)~=LOSindex
        counterMP(2)=counterMP(2)+1;
    end
end
disp(['Ho contato ' num2str(counterMP(2)) ' MP points su ' num2str(durationMP) ' totali (noise+MP).'])

%%% se il file è esteso a 40s di dataset
if sogliaMP1~=sogliaMP2
    
    for index=sogliaMP1+durationMP:sogliaMP2-1
        if indices(index)~=LOSindex
            counterMP(3)=counterMP(3)+1;
        end
    end
    disp(['Ho contato ' num2str(counterMP(2)) ' MP points su ' num2str(sogliaMP2-1-(sogliaMP1+durationMP)+1) ' totali (only noise).'])
    
    %%%
    for index=sogliaMP2:sogliaMP2+durationMP-1
        if indices(index)~=LOSindex
            counterMP(4)=counterMP(4)+1;
        end
    end
    disp(['Ho contato ' num2str(counterMP(4)) ' MP points su ' num2str(durationMP) ' totali (noise+MP).'])
    
    %%% Ultimi campioni di solo rumore
    for index=sogliaMP2+durationMP:Npoints
        if indices(index)~=LOSindex
            counterMP(5)=counterMP(5)+1;
        end
    end
    disp(['Ho contato ' num2str(counterMP(5)) ' MP points su ' num2str(Npoints-(sogliaMP2+durationMP)+1) ' totali (only noise).'])   
    
    %%%
    soglie=[sogliaMP1 sogliaMP1+durationMP-1 sogliaMP2 Npoints];
    results=struct('Npoints',Npoints,'counterMP',counterMP,'soglie',soglie,'durationMP',durationMP);
    
else
    soglie=[sogliaMP1 sogliaMP1+durationMP-1];
    results=struct('Npoints',Npoints,'counterMP',counterMP(1:2),'BeginMP',soglie,'durationMP',durationMP);
    
end
disp('-----------------------------------------------------------------------------------------------')

end