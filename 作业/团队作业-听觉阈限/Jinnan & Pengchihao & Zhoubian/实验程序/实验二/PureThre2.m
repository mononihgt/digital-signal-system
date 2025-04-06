%% Pure-tuned Threshold by JN 2024/04/08 

function PureThre2(subjectid)
sca;
try
    if nargin<1                              
        %如果没有输入参数的话
        subjectid='test';
        
    end
    %% screem normalization 屏幕校正；设置优先级为1：high priority level；不显示警告，隐藏鼠标
    %%%%personal computer no need.
    Priority(1);                                     % high priority level
    warning('off','all');                            % No warning
    %HideCursor;       
    WaitSecs(2);
    

    %% 设定屏幕
    background=128; %gray level of background 
    Screen('Preference','SkipSyncTests',1);
    [window,windowRect]=Screen('OpenWindow', 0, background); %,[0 0 500 500]
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    %% 设定阶梯法参数
    staircasetype='1u1d'; 
    strengthmatrix=[];
    responsematrix=[]; 
    minstrength=0.04;                                     
    maxstrength=1;
    initialstrength(1,:)=[minstrength(1) maxstrength(1)];
  
    %% 按键设置
    KbName('UnifyKeyNames');
    yeskey=KbName('j');
    nokey=KbName('k');
    qkey=KbName('q');
    %nexttrial=KbName('space');


    %% Basic Definition
    Fs = 5e3; % 采样率
    Duration_Pure = 0.100;          % 信号持续时间（秒）
    Duration_Noise = 1.000;         % 噪声持续时间（秒）
    t_Pure = 0:1/Fs:(Duration_Pure-1/Fs); % 纯音时间点
    Fre_Pure=1e3;                   % 纯音的频率
    half_band = [300 400 500];

    %% 实验开始
    for kk=[0 1]
        if kk==1        %% formal exp
            [paramatrix,~]=buildmatrix2;
            while KbCheck;  end                             % 清空键盘缓冲区,等待键盘上没有任何键被按下。
            Screen('DrawText',window,'Press any key to begin the formal exp',...
                round(windowRect(3)/2)-150,round(windowRect(4)/2),[0 0 0]);%提示开始   
            Screen('Flip', window);   
            KbWait;
            while KbCheck;  end  
        else            %% exercise
            [~,paramatrix]=buildmatrix2;
            while KbCheck;  end                             % 清空键盘缓冲区,等待键盘上没有任何键被按下。
            %Screen('FillOval',window,FixLum,FixationRect);  % 画注视点
            Screen('DrawText',window,'Press any key to begin some exercises',...
                round(windowRect(3)/2)-150,round(windowRect(4)/2),[0 0 0]);%提示开始   
            Screen('Flip', window);   
            KbWait;
            while KbCheck;  end  
        end
        
        for k=1:length(paramatrix(:,1))
           %% 每80试次休息一下
            if (mod(k,80)==1 && k~=1)
                 Screen('DrawText',window,'Please have a rest...Press any key to continue', ...
                     round(windowRect(3)/2)-150,round(windowRect(4)/2),[0 0 0]);
                 Screen('Flip', window);
                 WaitSecs(2);
                 while KbCheck; end
                 KbWait;
                 while KbCheck;  end                             % 清空键盘缓冲区,等待键盘上没有任何键被按下。          
            end
    
            %% 确定当前刺激强度(纯音/噪音)
            %%%%%%%%%   提取当前对应的刺激强度矩阵和反应矩阵
            Noise_index=paramatrix(k,2);     % the index of the noise   
            staircaseNumber=paramatrix(k,6);    % staircase number
            %%%%%%%%%   找到当前试次之前相同condition（noise x staircase）
            if k~=1
               strengthmatrix=paramatrix(find((paramatrix(1:k-1,2)==Noise_index)&(paramatrix(1:k-1,6)==staircaseNumber)),5);
               responsematrix=paramatrix(find((paramatrix(1:k-1,2)==Noise_index)&(paramatrix(1:k-1,6)==staircaseNumber)),7);
            end
            %%%%%%%%%   阶梯法,确定当前试次的刺激强度（呈现时间）
            PN_Ratio=updownstaircase(staircasetype,strengthmatrix,responsematrix,...
                    minstrength,maxstrength,initialstrength(mod(staircaseNumber,2)+1));
            
            if kk==0
                PN_Ratio=maxstrength;
            end

            %% Define noise
            
            Noise_filter=half_band(1,Noise_index);
            
            % % wgn函数需要Matlab中安装有Simulink APP                    
            % power_noise_dBW=10*log10(Noise_power);  % 转化为dBW单位
            % white_noise = wgn(1,Fs*Duration_Noise, power_noise_dBW);   
            
            % 若无法使用wgn函数，则使用randn函数。
            white_noise =sqrt(1)*randn(1, Fs*Duration_Noise);
            white_noise=filter_bandpass(Fs,white_noise,Noise_filter);
            
            Noise_power=mean(white_noise.^2);
            %% Define Pure
            Pure_power=PN_Ratio*Noise_power;
            A=sqrt(2*Pure_power);          % A为电压幅值，则正弦信号的功率为A^2/2
            pure_tone = A*sin(2*pi*Fre_Pure*t_Pure);
            
    
             %%  将纯音和白噪音叠加
            Stimuli_onsest=0.2*rand()+0.4;  %刺激在0.4~0.60s之中某一时刻出现
            Mute_pre=zeros(1, round(Fs*Stimuli_onsest));
            Mute_post=zeros(1,length(white_noise)-length(Mute_pre)-length(pure_tone));
            Time_Align_signal = [Mute_pre, pure_tone, Mute_post];
            combined_signal =Time_Align_signal +white_noise;
            
            combined_signal=[combined_signal;combined_signal];%% 双声道
            
            %% 试次开始
            Screen('DrawText',window,'+',round(windowRect(3)/2),round(windowRect(4)/2),[0 0 0]);
            Screen('DrawText',window,'Plz listen !',round(windowRect(3)/2)-50,round(windowRect(4)/3),[0 0 0]);
            Screen('Flip', window);
            WaitSecs(0.5);
            
         
            %% 播放合并后的信号
            sound(combined_signal, Fs);
            WaitSecs(length(combined_signal(1,:))/Fs+0.1);
    
            %% 呈现刺激后注视点
            Screen('DrawText',window,'+',round(windowRect(3)/2),round(windowRect(4)/2),[0 0 0]);
            Screen('DrawText',window,'Did you hear a pure tone?',round(windowRect(3)/2)-120,round(windowRect(4)/4),[0 0 0]);
            Screen('DrawText',window,'[J] for Yes; [K] for No',round(windowRect(3)/2)-110,round(windowRect(4)/4)+40,[0 0 0]);
            secsonset=Screen('Flip', window);
    
            %% 按键
            keycode=zeros(1,256);
            while (keycode(qkey)==0)&&(keycode(yeskey)==0)&&(keycode(nokey)==0)
                WaitSecs(0.0001);
                [~,secs,keycode]=KbCheck;
            end
            RT=secs-secsonset;
            if keycode(qkey)==1 %q键退出
                break
            end  
            
            %% 记录反应
            
            if (keycode(yeskey)==1 && sum(keycode)==1)         %%k按键表示看到向上运动
                Response=1;
            else
                Response=0;
            end
            
            paramatrix(k,3)=Noise_power;
            paramatrix(k,4)=Pure_power;
            paramatrix(k,5)=PN_Ratio;
            paramatrix(k,7)=Response;
            paramatrix(k,8)=RT;
            paramatrix(k,9)=Stimuli_onsest;
            
            if  kk==0
                save([subjectid '_practice_paramatrix'],"paramatrix");
            else
                save([subjectid '_formal_paramatrix'],"paramatrix");
            end
            
            %% 等待进入下一个试次
            Screen('DrawText',window,'+',round(windowRect(3)/2),round(windowRect(4)/2),[1 1 1]);
            Screen('Flip', window);
            WaitSecs(0.4+0.2*rand());
            while KbCheck;  end 
        end
        if keycode(qkey)==1 %q键退出
                break
        end 
    end
    Screen('DrawText',window,'The exp is over! Plz wait for 10 seconds for you results',round(windowRect(3)/2)-200,round(windowRect(4)/4),[0 0 0]);
    Screen('Flip', window);
    fittingresult(subjectid,3,[]);
    ShowCursor;
    sca;

catch
    ShowCursor;
    sca;
    psychrethrow(psychlasterror);
end

end
