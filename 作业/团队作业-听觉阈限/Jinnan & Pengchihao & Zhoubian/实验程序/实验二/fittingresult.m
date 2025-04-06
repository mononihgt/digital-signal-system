%   横轴为刺激强度
%   纵轴为判断正确概率
%   最后图形为每个刺激强度下，正确判断概率的拟合曲线图
%   标出阶梯法对应概率时刺激强度。1u1d-50%;  1u2d-70.71%;  1u3d-79.37%
%   处理数据

function threshold=fittingresult(subjectid,condition_num,dotdeleted)
    subjectid_temp=subjectid;
    paramatrix=[];
    threshold=[];
    subjectid=[subjectid_temp '_formal'];
    load([subjectid '_paramatrix']);
    
    rawdata=[];
    StrthRsp=[];
    for k=1:condition_num
        rawdata(:,1)=paramatrix(paramatrix(:,2)==k ,5); %% stim strength
        rawdata(:,2)=paramatrix(paramatrix(:,2)==k ,7); %% response
        rawdata=round(rawdata*100000)/100000;
        StrthRsp(:,1)= unique(rawdata(:,1)); %唯一化
        for i=1:length(StrthRsp(:,1))
            StrthRsp(i,2)=length(find(rawdata(:,1)==StrthRsp(i,1)&rawdata(:,2)==1));%每个刺激强度的正确判断的试次数(block)
            StrthRsp(i,3)=length(find(rawdata(:,1)==StrthRsp(i,1)));                %每个刺激强度的总的试次数
        end
        StrthRsp(:,4)=StrthRsp(:,2)./StrthRsp(:,3);
        StrthRsp(dotdeleted,:)=[];
        StimStrength=StrthRsp(:,1);     %唯一化的各刺激强度
        Respright=StrthRsp(:,4);        %判断的正确率
        TrialCount=StrthRsp(:,3);       %每个刺激条件总的试次数
        threshold(1,k)=fittingsigmoid(StimStrength,Respright,TrialCount); 
        figure_name=[subjectid '_' num2str(k)];
        saveas(gcf,[figure_name,'.jpg']);
        close all
        rawdata=[];
        StrthRsp=[];
    end
    save([subjectid ' Pure-tone Threshold'],"threshold")
end



%%该函数用于s形曲线拟合
%%StimStrength α值
%%Respright     正确判断的概率
%%TrialCount    每个刺激点的trial数

function Threshold=fittingsigmoid(StimStrength,Respright,TrialCount)
    %function output=fittingsigmoid(StimStrength,Respright,TrialCount,subjectid,condition)
    
    
    %%%%%%%%%%%%%%拟合高斯分布曲线%%%%%%%%%%%%%%%%
    x = log10(min(StimStrength)):0.001:log10(max(StimStrength));
    %Threshold0= (log10(min(StrthRsp(:,1)))+ log10(max(StrthRsp(:,1))))/2;
    tempindex=find(TrialCount==max(TrialCount));
    Threshold0=log10(StimStrength(tempindex(1)));%以trial数最多的强度作为拟合起始参数
    %Threshold0=log10(4);%以trial数最多的强度作为拟合起始参数
    slope0=0.3;
    %amp0=3;
    beta0 = [Threshold0 slope0];
    fun = @(beta,x)(1+erf(((x-beta(1))/beta(2))/sqrt(2)))/2;                %%%chance level=0    情况
    %fun = @(beta,x)((1+erf(((x-beta(1))/beta(2))/sqrt(2)))/2)*0.5+0.5;     %%%chance level=0.5  情况
    %fun = @(beta,x)((1+erf(((x-beta(1))/beta(2))/sqrt(2)))/2)*0.75+0.25;   %%%chance level=0.25 情况
    
    %2014版本matlab拟合命令
    wnlm = fitnlm(log10(StimStrength),Respright,fun,beta0,'weight',TrialCount);
    r_square=wnlm.Rsquared.Ordinary;
    yfit = predict(wnlm,x');
    % SSR=sum(yfit-mean(Respright))*sum(yfit-mean(Respright));
    % SST=sum(Respright-mean(Respright))*sum(Respright-mean(Respright));
    % R2=SSR/SST;
    %%2011或之前版本matlab拟合命令
    % [wnlm r2] = nlinfit(log10(StimStrength),Respright,fun,beta0);
    % StimStrength
    % RespAccuracy
    % TrialCount
    % yfit=fun(wnlm,x');
    
    %%%%%%%%%%%%%%画 图%%%%%%%%%%%%%%%%%%%%%%%%%
    %图形的颜色
    
    linecolorr=[0 0 1];
    dotcolorr=[0.5 0.5 0.8];
    
    
    %画出点、线
    hold on
    plot(10.^x,yfit,'LineWidth',2,'Color',linecolorr)
    scatter(StimStrength,Respright,TrialCount*5,'ro','MarkerFaceColor',dotcolorr,'MarkerEdgeColor','none')
    ThreP=0.5;%1u1d
    ThreIndex=abs(yfit-ThreP)==min(abs(yfit-ThreP));
    Threshold=10^x(ThreIndex);%找到对应的阈值
    line([min(StimStrength) Threshold],[ThreP ThreP],'color',[0.5 0.5 0.5])
    line([Threshold Threshold],[ThreP 0],'color',dotcolorr)
    text(Threshold,0.05,num2str(Threshold,'%4.3g'));
    xlabel('Power_{Pure tone}/Power_{Noise} (dB)')
    ylabel('Probability')

end
