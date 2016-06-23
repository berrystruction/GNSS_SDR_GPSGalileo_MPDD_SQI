%function [ detectprofile, judmentfunc ] = MPDD_part2(MPDD_part1,CN0profile,Tepoch,mvavgtime,minCN0considered,settings)
function [ detectprofile ] = MPDD_part2(MPDD_part1,CN0profile,settings,settingsMPDD,timesize,DistortionLAFindex)

Tepoch=settings.navSolPeriod;
mvavgtime=settingsMPDD.mvgAvgtime;
minCN0considered=settingsMPDD.minCN0considered;
%MPDD_PART2 Summary of this function goes here
% minCN0considered=34;
lenepoch=floor(Tepoch/mvavgtime);
%epochDuration=(1:lenepoch);
epochDuration=linspace(0,Tepoch*1e-3,lenepoch);

%seriousness=normpdf(epochDuration,lenepoch/2,lenepoch/2/10).^2;
seriousness=normpdf(epochDuration,1,lenepoch/2/60).^2;
seriousness=seriousness/max(seriousness)+1; %seriousness=ones(1,Navg);
%half=lenepoch/2; seriousness=[exp(epochDuration(1:half)) exp(-epochDuration(half+1:end))];
%xbins=[-1 1];
MPDDoffset=settingsMPDD.MPDDoffset; % [MPDD point] 10 is about 500 ms (offset to put istant of PVT publication in the middle of the MPDD epoch)

%% Section 2 - Detection
index1=MPDDoffset; index2=index1+lenepoch-1; index=1;
detectprofile=zeros(1,floor(settings.msToProcess/1e3));


% Two kind of loop depending what type of LAF used
if nargin==6
    % Type 4 - settingsMPDD.LAFversion=4 - LS
    while index2<timesize
        pieceOfCN0=CN0profile(index1:index2);
        %meanCN0=mean(pieceOfCN0);
        CN0function=(1+(atan(pieceOfCN0-minCN0considered)/(pi/2)))/2;
        %size(seriousness), size(CN0function), size(MPDD_part1(index1:index2))
        judmentfunc=seriousness.*CN0function.*MPDD_part1(index1:index2).*DistortionLAFindex(index1:index2);
        detectprofile(index)=sum(judmentfunc);
        
        index1=index2+1;
        index2=index2+lenepoch;
        index=index+1;
    end
    
else
    % Type 3 - settingsMPDD.LAFversion=3 - CLS and others
    while index2<timesize+1
        pieceOfCN0=CN0profile(index1:index2);
        CN0function=(1+(atan(pieceOfCN0-minCN0considered)/(pi/2)))/2;
        
        judmentfunc=seriousness.*CN0function.*MPDD_part1(index1:index2);
        detectprofile(index)=sum(judmentfunc);
        
        index1=index2+1;
        index2=index2+lenepoch;
        index=index+1;
    end
    
end





% normalization between -1 and 1
maxx=sum(seriousness*ones(settingsMPDD.MPDDNpoints,1)); % *1 all measurement are made under high CN0 hypothesis
detectprofile=(detectprofile/maxx+1)/2;

% index1=(index-1)*lenepoch+MPDDoffset;
% index2=index1+lenepoch-1;

% figure
% subplot(221)
% plot(epochDuration,judmentfunc,'.-',epochDuration,sum(judmentfunc)*ones(1,lenepoch),'-',epochDuration,mean(judmentfunc)*ones(1,lenepoch),'-'), legend('judment function','Sum','Average')
% title('Final judment function')
% grid on
% subplot(222)
% plot(epochDuration,MPDD_part1(index1:index2),'o')
% title('MPDD results')
% grid on
% subplot(223)
% [g, xg]=hist(judmentfunc,2);
% bar(xg,g/Navg*100)
% title('Final Decision')
% grid on
% subplot(224)
% [f, xf]=hist(MPDD_part1,2);
% bar(xbins,f/Navg*100)
% title('Final Decision without judment function')
% grid on


end

