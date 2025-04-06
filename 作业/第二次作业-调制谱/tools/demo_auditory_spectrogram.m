% This code illustrates how to compute and plot the
% auditory spectrogram simulated by a cochlear model
% reference: Yang, Wang, & Shamma. IEEE Trans Information Theory 38.2 (1992): 824-839.
% author: Nai Ding, ding_nai@zju.edu.cn

clear;clc;
% load the sound file, resample it to 16 kHz if necessary
[x,fs]=audioread('/Users/maopeixuan/Documents/大三/大三下/信号与认知系统/第二次作业/materials/dialect.mp3');
x=resample(x,16e3,fs);
% generate the auditory spectrogram
[v,cf]=wav2aud2(x,[5 8 -2 0]);
% plot the auditory spectrogram
figure;imagesc([1 size(v,1)]/200,1:128,20*log10(max(v,0))');axis xy
set(gca,'ytick',1:30:128);set(gca,'yticklabel',cf(1:30:128));
xlabel('time (s)');ylabel('frequency (Hz)');title('auditory spectrogram')
l=colorbar;title(l,'dB');
save('materials/dialect.mat','v');
save('materials/dialect_cf.mat','cf');

