
clc;clear;

%根据设置频率与持续时间的参数
                
%melody=[1 4 5 6 7 6 5 6 5 4 5 2 4 4 5 3 2 1 2 ]  简谱
melodyHalfTone=[-2 3 5 7 9 7 5 7 5 3 5 0 3 3 5 2 0 -2 0];% 简谱音符距离D4的相对半音数
duration0=[0.25 0.5 0.125 0.125 0.5 0.125 0.125 0.5 0.125 0.125 0.25 0.5 0.5 0.125 0.125 0.25 0.375 0.125 0.3];%每个音持续时间
duration=duration0*2.5;
frequence=440*2.^((melodyHalfTone)/12);% 每个音符的频率

% 设置采样率和总持续时间
fs = 44100; % 采样率（每秒的样本数）


% 生成音乐

music = [];
for note = 1:length(frequence)
    duration_note = duration(note);  
    t = 0:1/fs:duration_note;

    wave = 0.8*sin(2*pi*frequence(note)*t); % 初始化波形

   

    %动态调整谐波

    if frequence(note) < 293.66 % D4以下
        dr=2.8;
        for h = 2:5 
            amplitude = (1 / h^dr); 
            harmonic = amplitude * sin(2*pi*frequence(note)*h*t);
            wave = wave + harmonic;
        end
    elseif frequence(note) >= 293.66 && frequence(note) < 587.33 % D4到D5之间
        dr=2.5;
        for h = 2:6 
            amplitude = (1 / h^dr); 
            harmonic = amplitude * sin(2*pi*frequence(note)*h*t);
            wave = wave + harmonic;
        end
    else
        dr=2.2;
        for h = 2:7
            amplitude = (1 / h^dr); 
            harmonic = amplitude * sin(2*pi*frequence(note)*h*t);
            wave = wave + harmonic;
        end
    end

    
    wave = wave / max(abs(wave)); % 归一化

    % 动态共振滤波
    wave = dynamicResonanceFilter(wave, frequence(note)); 
     

%   动态调整ASDR包络
    ADSR = dynamicASDR(frequence(note));
    A = ADSR(1); % 攻击时间
    D = ADSR(2); % 衰减时间
    S = ADSR(3); % 维持水平
    R = ADSR(4)*duration_note*1.75; % 释放时间

    release_curve = linspace(S, 0, floor(R * fs)).^1.5; % 维持之前的衰减曲线

    % 构建完整的ADSR包络
    envelope = [linspace(0, 1, floor(A * fs )),...
                linspace(1, S, floor(D * fs)), ...
                S * ones(1, length(t) - length([linspace(0, 1, floor(A * fs)), linspace(1, S, floor(D * fs)), release_curve])), ...
                release_curve];
    envelope = envelope(1:length(t)); % 确保包络与音符长度一致
    wave = wave .* envelope;
    

    % 共振峰处理和失真处理
    dynamicResonance = @(noteFreq) [noteFreq * 0.9, noteFreq * 1.1; noteFreq * 1.8, noteFreq * 2.2; noteFreq * 2.7, noteFreq * 3.3];
    resonantFrequencies = dynamicResonance(frequence(note));
    filteredWave = wave;
        for i = 1:size(resonantFrequencies, 1)
            fc = mean(resonantFrequencies(i, :));
            bandwidth = diff(resonantFrequencies(i, :));
            filteredWave = filteredWave + 0.2 * bandpass(wave, [fc-bandwidth/2, fc+bandwidth/2], fs);
        end

    
    waveDistorted = tanh(filteredWave*0.7);

    music = [music waveDistorted];

end

%滤波以去除每个音最后的爆鸣声
fc = 100; % 截止频率（Hz）
Wn = fc/(fs/2); % 归一化截止频率
order = 8; % 滤波器阶数
b = fir1(order, Wn, 'low'); % 设计8阶低通FIR滤波器系数

% 应用滤波器
filtered_signal = filter(b, 1, music);

% 播放音乐
sound(music, fs);

windowlen = 5;
ytime = smoothdata(music,'movmean',windowlen);
yfrequency = fft(ytime,fs);
%%
yfrequency=abs(fftshift(fft(ytime)));           %FFT变换
figure(2); 
subplot(2,1,1);plot(ytime); title('滤波器时域图');       %时域波形
subplot(2,1,2);plot(yfrequency); title('滤波器频谱图 ');    %频谱图




% 动态ASDR包络函数
function ADSR = dynamicASDR(frequency)
    % 根据频率动态调整ADSR参数
    if frequency < 293.66 % D4以下
        A = 0.013;
        D = 0.1;
        S = 0.5;
        R = 0.6;
    elseif frequency >= 293.66 && frequency < 587.33 % D4到D5之间
        A = 0.01;
        D = 0.1;
        S = 0.5;
        R = 0.55;
    else % D5以上
        A = 0.007;
        D = 0.1;
        S = 0.5;
        R = 0.5;
    end
    ADSR = [A, D, S, R];

end


% 动态共振滤波器函数
function filteredWave = dynamicResonanceFilter(wave, freq)
    if freq < 500 % 低频区域
        b = [0.1, 0.2, 0.1]; % 低频强化
    else % 高频区域
        b = [0.06, 0.1, 0.06]; % 高频强化
    end
    a = [1, -0.4, 0.2]; % 维持共振体特性
    filteredWave = filter(b, a, wave); % 应用滤波器
end
