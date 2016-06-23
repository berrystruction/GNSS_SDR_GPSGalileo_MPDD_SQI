function [nbits]=test_raw_samples(fbin,sampfreq,seek_sec)


% clear all;
% close all;         
% % sampfreq = 39.8e6;
% sampfreq = 16.3676e6; % sampling frequency [Hz]
% 
% % Test Front end
% 
% % open the input file with raw samples
% % fbin = fopen('E:\GNSS SIGNALS - DATA SETS\Agilent Hardware Generator\santarosa.bin','rb');
% %fbin = fopen('.\prova.bin','rb');
% %fbin = fopen('.\ssdm2_GPSPRN1_pulse500_off10.bin','rb');
% 
% % Select the log file
% [filename, pathname] = uigetfile({'*.bin;*.grab'}, 'Select a file');
% name = fullfile(pathname, filename)
% 
% fbin = fopen(name,'rb');

%global seek_sec;

% read a snapshot of sample with appropriated format
status = fseek(fbin,ceil(sampfreq*seek_sec),-1);
% status = fseek(fbin,ceil(sampfreq*0.1),-1);
%status = fseek(fbin,ceil(sampfreq*15.1),-1); % Skip 15 s

if (status)
    error('fseek error')
end
 
% [data, cnt]= fread(fbin,10e6,'uchar'); % Campionato con degli uchar tra 0 e 255!!
[data, cnt]= fread(fbin,10e6,'int8'); % Campionato con degli schar 
mean_value = mean(data)
% data = data-127; % Campionato con degli uchar tra 0 e 255!!

% fclose(fbin);

% plot the first 1000 samples
figure;
plot(data(1:1000),'bo-'),grid on;box on;
axis tight;aa=axis;axis([aa(1) aa(2) aa(3)-1 aa(4)+4])
title('Time Domain of the first 1000 samples');
xlabel('Num. of samples');ylabel('Amplitude');

x = linspace(-128,127,256);
figure,counts=hist(data(1:1e6),x);
bar(x,counts)
grid on;box on;

title('Histogram of received samples')
% tento di capire il numero di bit cercando i min valore delle x con un
% valore diverso da 0
indexmin=0;
flagMinx=false;
while(flagMinx==false)
    indexmin=indexmin+1;
    if counts(indexmin)~=0
        flagMinx=true;
    end
end
indexmax=length(x);
flagMax=false;
while(flagMax==false)
    indexmax=indexmax-1;
    if counts(indexmax)~=0
        flagMax=true;
    end
end
nbits=ceil(log2(abs(max(abs(x(indexmin)),abs(x(indexmax))))))+1;


% plot signal spectrum
% numptsfft = 1048576;  
numptsfft = length(data);    
NFFT = 4096*4; 

[Pxx_L1,F_L1] = psd(data(100:numptsfft)-mean(data(100:numptsfft)),NFFT,sampfreq);
figure
plot(F_L1/1e6,10*log10(abs(Pxx_L1)),'b'),grid on;box on;
axis tight;aa=axis;axis([aa(1) aa(2) aa(3)-1 aa(4)+4])
% title('L1 Spectrum');
title('Spectrum');
xlabel('Frequeny [MHz]');ylabel('Magnitude [dB]');

end
