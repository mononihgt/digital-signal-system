clear; clc; close all;

% --- 添加路径并加载数据 ---
addpath ./tools % 假设您的bootstrap函数在这个文件夹

% 加载EEG数据
load Data/ybE1 %% EEG data in Experiment 1 (assuming 'yb' is the variable name)
% 或者如果您保存的文件名不同，用下面这行:
% load Data/脑电数据 %% EEG data in Experiment 1 (assuming 'x' is the variable name)
% yb = x; % 如果变量名是x

% --- 参数设置 ---
fs = 16; % 降采样后的频率 (Hz)
t = 0:(1/fs):(size(yb,1)/fs - 1/fs); % 创建时间向量
t = t - t(find(t>=1, 1)-1); % 调整时间，使 t=1 对应分析窗开始

wind = (t >= 1 & t < 11); % 分析窗口：1到11秒 (共10秒)
TD = sum(wind); % 分析窗口内的采样点数
channel_idx = 48; % Cz通道的索引 (请确认!)
[~, ~, num_trials, num_subjects] = size(yb); % 获取维度信息

f = (0:TD-1)' / TD * fs; % 计算频率轴
freq_range_plot = (f >= 0.5 & f <= 5); % 绘图时关注的频率范围
num_repetitions_subj = 10; % 被试抽样的重复次数（用于第二部分分析）
num_bootstrap_trial = 100; % 每个试次内进行bootstrap的次数

rng('default'); % 为了结果可复现，固定随机数种子

% --- 辅助函数：计算功率谱 ---
% 输入：data_in (时间点, 实体数) - 实体可以是试次或被试
% 输出：power_spectrum (频率点, 实体数)
function power_spec = calculate_power_spectrum(data_in, TD, fs)
    fft_result = fft(data_in, TD);
    power_spec = (2 * abs(fft_result(1:TD,:)) / TD).^2;
end

% --- 1. 分析试次数目对频谱的影响 (每个试次单独分析, 使用Bootstrap) ---
fprintf('Analyzing effect of individual trials using Bootstrap...\n');
power_vs_trials = zeros(TD, num_trials); % 预分配最终结果矩阵 (频率点 x 试次数)

% 选择Cz通道数据
cz_data_all_subjects = squeeze(yb(wind, channel_idx, :, :)); % (时间点, 试次, 被试)

for k = 1:num_trials
    fprintf('  Processing trial %d:\n', k);
    % 选择第 k 个试次，所有被试的数据
    data_trial_k = squeeze(cz_data_all_subjects(:, k, :)); % (时间点, 被试)
    if num_subjects == 1 % 处理单被试情况
       data_trial_k = reshape(data_trial_k, TD, 1);
    end

    temp_spectra_trial_k = zeros(TD, num_bootstrap_trial); % 存储当前试次的Bootstrap结果

    for b = 1:num_bootstrap_trial
        % 1. 有放回地抽取 num_subjects 个被试索引
        resample_indices = randi(num_subjects, num_subjects, 1);

        % 2. 提取对应被试的数据
        resampled_data = data_trial_k(:, resample_indices);

        % 3. 计算每个重采样被试的功率谱
        power_per_resampled_subject = calculate_power_spectrum(resampled_data, TD, fs);

        % 4. 计算这次bootstrap样本的几何平均功率谱
        bootstrap_geom_mean = geomean(power_per_resampled_subject, 2); % 跨重采样被试几何平均

        % 5. 存储本次bootstrap的结果
        temp_spectra_trial_k(:, b) = bootstrap_geom_mean;

        % (可选) 显示Bootstrap进度
        % if mod(b, 20) == 0
        %     fprintf('    Bootstrap rep %d/%d\n', b, num_bootstrap_trial);
        % end
    end

    % 计算 num_bootstrap_trial 次bootstrap结果的算术平均值，作为试次k的最终频谱
    power_vs_trials(:, k) = mean(temp_spectra_trial_k, 2);

    fprintf('  Trial %d done.\n', k);
end

% --- 绘制试次数影响的图像 (Bootstrap估计) ---
figure;
imagesc(1:num_trials, f(freq_range_plot), 10*log10(power_vs_trials(freq_range_plot, :))); % dB单位
set(gca, 'YDir', 'normal');
xlabel('Trial Number');
ylabel('Frequency (Hz)');
title(sprintf('EEG Spectrum (Cz) per Trial (Bootstrap Avg across %d Subjects, %d Reps)', num_subjects, num_bootstrap_trial));
colorbar;
colormap('jet');


% --- 2. 分析被试数量对频谱的影响 (随机抽样+重复) ---
% (这部分代码保持不变)
fprintf('\nAnalyzing effect of number of subjects (random subsampling)...\n');
power_vs_subjects = zeros(TD, num_subjects); % 预分配最终结果矩阵

% 先计算好每个被试的试次平均信号
cz_data_subject_avg = squeeze(mean(yb(wind, channel_idx, :, :), 3)); % (时间点, 被试)

for m = 1:num_subjects
    fprintf('  Processing for %d subjects:\n', m);
    temp_spectra_for_m = zeros(TD, num_repetitions_subj); % 存储当前m值下的重复结果

    for rep = 1:num_repetitions_subj
        % 随机抽取 m 个被试的索引（无放回）
        selected_subject_indices = randperm(num_subjects, m);

        % 提取这 m 个被试的试次平均数据
        data_m_selected = cz_data_subject_avg(:, selected_subject_indices);

        % 计算这 m 个被试各自的功率谱
        power_per_subject_m = calculate_power_spectrum(data_m_selected, TD, fs);

        % 计算这 m 个被试的几何平均功率谱
        group_power_spectrum_rep = geomean(power_per_subject_m, 2); % 跨选定被试几何平均

        temp_spectra_for_m(:, rep) = group_power_spectrum_rep; % 存储本次重复的结果

        fprintf('    Repetition %d done.\n', rep);
    end

    % 计算重复次数的算术平均值，作为最终m个被试的结果
    power_vs_subjects(:, m) = mean(temp_spectra_for_m, 2);

end

% --- 绘制被试数量影响的图像 (随机抽样+重复) ---
figure;
imagesc(1:num_subjects, f(freq_range_plot), 10*log10(power_vs_subjects(freq_range_plot, :))); % dB单位
set(gca, 'YDir', 'normal');
xlabel('Number of Subjects Included (Randomly Sampled)');
ylabel('Frequency (Hz)');
title(sprintf('Effect of Number of Subjects on Group EEG Spectrum (Cz, %d Reps Avg)', num_repetitions_subj));
colorbar;
colormap('jet');

fprintf('\nAnalysis complete.\n');