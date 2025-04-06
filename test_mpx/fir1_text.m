clc;clear;close all;
% 参数设置
fs = 16000; % 采样频率，单位为 Hz
A = 0.2; % 示例值，根据实际需求调整
f1 = (1-A)*1000; % 带通滤波器下限频率
f2 = (1+A)*1000; % 带通滤波器上限频率
n = 50; % 滤波器阶数，可根据需求调整

% 归一化频率
Wn = [f1 f2] / (fs/2);

% 设计 FIR 带通滤波器
b = fir1(n, Wn, 'bandpass');

% 生成白噪声
t = 0:1/fs:1; % 时间序列
noise = randn(size(t)); % 生成白噪声

% 应用滤波器
filtered_noise = filter(b, 1, noise);

% FFT 参数
N = length(t); % FFT 点数
f = (0:N-1)*(fs/N); % 频率轴

% 计算 FFT
fft_noise = fft(noise);
fft_filtered_noise = fft(filtered_noise);

% 计算幅度谱（取模并归一化）
fft_noise_mag = abs(fft_noise) / N;
fft_filtered_noise_mag = abs(fft_filtered_noise) / N;

% 绘制频谱图
figure;
subplot(2,1,1);
plot(f, fft_noise_mag);
title('原始白噪声的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
xlim([0 fs/2]); % 只显示正频率部分

subplot(2,1,2);
plot(f, fft_filtered_noise_mag);
title('滤波后信号的频谱');
xlabel('频率 (Hz)');
ylabel('幅度');
xlim([0 fs/2]); % 只显示正频率部分

% 使用 freqz 分析滤波器的频率和相位响应
figure;
freqz(b, 1, 1024, fs); % 绘制频率响应和相位响应

% 如果需要单独绘制幅度响应和相位响应
[h, f] = freqz(b, 1, 1024, fs);

figure;
subplot(2,1,1);
stem(f, 20*log10(abs(h))); % 幅度响应（分贝）
title('滤波器的幅度响应');
xlabel('频率 (Hz)');
ylabel('幅度 (dB)');

subplot(2,1,2);
stem(f, unwrap(angle(h))*180/pi); % 相位响应（度）
title('滤波器的相位响应');
xlabel('频率 (Hz)');
ylabel('相位 (度)');