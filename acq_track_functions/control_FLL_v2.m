
% Control FLL.
% This program control the FLL status, checking the freqeuncy of the local
% carrier frequency, passed from the FLL program.
% The check is based on a FIFO structure and is performed at a predefined
% control time, which could be for example every 0.25 s.

function [FLL_status,freq_control,Freq_sum_out,Old_Freq_sum_out] = control_FLL_v2(act_freq,loop,thres,Freq_sum,Old_Freq_sum,T_int)


% global Freq_sum;
% global Old_Freq_sum;
% global T_int

Freq_sum = (1/20)*act_freq + ((20-1)/20)*Freq_sum;

FLL_status = 0;
freq_control = 0;

if (rem(loop*(T_int*1e3),100)==0)
    
    if (abs((Freq_sum - Old_Freq_sum))< thres)
        FLL_status = 1;
        freq_control = Freq_sum;
    end
   % figure(100),hold on, plot(loop, Freq_sum,'o');
    Old_Freq_sum = Freq_sum;
end

Old_Freq_sum_out=Old_Freq_sum;
Freq_sum_out=Freq_sum;
% if(loop >= 5000)       % After 5 seconds (5000 msec, since the PDI used in the FLL is equal to 1 msec), the time is out! FLL not locked
% if(loop >= 10000)       % After 10 seconds (10000 msec, since the PDI used in the FLL is equal to 20 msec), the time is out! FLL not locked
%     FLL_status = -1;
% end

end