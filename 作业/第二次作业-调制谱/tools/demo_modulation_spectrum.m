% This code illustrate how to compute and plot the modulation spectrum.
% The input audio is segmented into 6-second chunks and the modulation
% spectrum is calculated for each chunk and averaged.
% reference: Ding, Patel, & Poeppel, to be published
% author: Nai Ding, gahding@gmail.com

clear;clc;close all
% specify the name of the audio file to analyze
filename='/Users/maopeixuan/Documents/大三/大三下/信号与认知系统/第二次作业/materials/cello.mp3';

% get the sampling rate, fs
ainfo=audioinfo(filename);
fs=ainfo.SampleRate;
% load 6-seconds duration non-overlapping chunks of the audio file
for ind=1:1000
  % load a 6-seconds duration chunk
  try  [x,fs]=audioread(filename,[fs*6*(ind-1)+1 fs*6*ind]);
  catch
    if exist('ms')
      break;
    else
      [x,fs]=audioread(filename);
      x(6*fs,1)=0;
    end;
  end
  if size(x,2)>1,x=mean(x,2);end
  disp(['...' ' processing chunk ' num2str(ind)])
  % resample the recordign to 16 k Hz
  x=resample(x,16e3,fs);
  % generate the auditory spectrogram
  [v,cf]=wav2aud2(x,[5 8 -2 0]);
  % compute the narrow-band modulation spectrum
  [ms(:,ind),f]=Modulation_Spectrum(v,200,'log');
end

  
disp(['The analyzed audio file is segmented into ' num2str(ind) ' chunks.'])

% calculate the root-mean-square of the modulation spectrum
% across all the 6-second duration chunks
ms_rms=sqrt(mean(ms.^2,2));
ms_rms=ms_rms/max(ms_rms([f<32 & f>=.5]));
% plot the modulation spectrum averaged over chunks
save('materials/cello_ms.mat','ms');
save('materials/cello_ms_rms.mat','ms_rms');
save('materials/cello_f.mat','f')
figure;
plot(f(1:size(ms,1)),ms_rms)
xlim([.5 32]);
set(gca,'xscale','log')
set(gca,'xtick',[.5 1 2 4 8 16 32])
xlabel('frequency (Hz)')
ylabel('normalized amplitude')