%% 20240408 Pure-tune Threshold by JN
%% column 1 trial number
%% COLUMN 2 噪声的功率序号
%% COLUMN 3 噪音的真实功率
%% COLUMN 4 纯音的真实功率
%% COLUMN 5 纯音的功率/噪音的功率
%% column 6 staircase的情况：1 2
%% column 7 判断情况：1，听得见纯音；0，听不见纯音
%% column 8 反应时

function [paramatrix,paramatrix_prac] = buildmatrix(noise_num)

power_noise=1:noise_num;        % 噪声功率水平的个数

StaircaseNum=[1 2];     % 每个condition有2个staircase
trialpercondition=40;   % 80个trials
trialpercondition_practice=1;
[x0,x1]=ndgrid(power_noise,StaircaseNum);   % 生成参数组合矩阵
combinedpara=[x0(:),x1(:)];                 % 实验trial所有的条件组合

col=8;  % 参数矩阵一共8列
paramatrixtemp0=zeros(length(combinedpara(:,1)),col);
paramatrixtemp0(:,2)=combinedpara(:,1);     % 噪声的功率
paramatrixtemp0(:,6)=combinedpara(:,2);     % 阶梯序号

%% fomal exp
paramatrix0=[];
for ii=1:trialpercondition
    paramatrix0=[paramatrix0;paramatrixtemp0];
end

paramatrix=[];

[pointer,index]=Shuffle(paramatrix0(:,1));%%%打乱刺激的次序

for i=1:length(paramatrix0(:,1))
    paramatrix(i,:)=paramatrix0(index(i),:);
end

paramatrix(:,1)=(1:length(paramatrix(:,1)))';

%% practice
paramatrix0_prac=[];
for ii=1:trialpercondition_practice
    paramatrix0_prac=[paramatrix0_prac;paramatrixtemp0];
end

paramatrix_prac=[];

[pointer,index]=Shuffle(paramatrix0_prac(:,1));%%%打乱刺激的次序

for i=1:length(paramatrix0_prac(:,1))
    paramatrix_prac(i,:)=paramatrix0_prac(index(i),:);
end

paramatrix_prac(:,1)=(1:length(paramatrix_prac(:,1)))';




end