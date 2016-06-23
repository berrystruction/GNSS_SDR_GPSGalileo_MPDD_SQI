function majorityIndex=majorityRule(vecIndex,vtime,PRN)
Nindex=size(vecIndex,1);

aus1=vecIndex>0;
%aus2=-1*(vecIndex<=1);
%newInd=aus1+aus2;
threshold=ceil(Nindex/2);
newInd=sum(aus1,1);
majorityIndex=newInd>=threshold;

majorityIndex=2*majorityIndex-1;
if 1==0
    figure
    subplot(211)
    plot(vtime,newInd,'o',vtime,ones(size(vtime))*threshold,'r-')
    grid on
    title(['majority rule for PRN ' num2str(PRN)])
    xlabel('Time')
    
    subplot(212)
    plot(vtime,majorityIndex,'o')
    title('Final results')
    xlabel('Time')
    grid on
end

end