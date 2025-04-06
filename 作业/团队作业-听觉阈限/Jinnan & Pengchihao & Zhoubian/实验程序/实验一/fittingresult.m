%   ����Ϊ�̼�ǿ��
%   ����Ϊ�ж���ȷ����
%   ���ͼ��Ϊÿ���̼�ǿ���£���ȷ�жϸ��ʵ��������ͼ
%   ������ݷ���Ӧ����ʱ�̼�ǿ�ȡ�1u1d-50%;  1u2d-70.71%;  1u3d-79.37%
%   ��������

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
        StrthRsp(:,1)= unique(rawdata(:,1)); %Ψһ��
        for i=1:length(StrthRsp(:,1))
            StrthRsp(i,2)=length(find(rawdata(:,1)==StrthRsp(i,1)&rawdata(:,2)==1));%ÿ���̼�ǿ�ȵ���ȷ�жϵ��Դ���(block)
            StrthRsp(i,3)=length(find(rawdata(:,1)==StrthRsp(i,1)));                %ÿ���̼�ǿ�ȵ��ܵ��Դ���
        end
        StrthRsp(:,4)=StrthRsp(:,2)./StrthRsp(:,3);
        StrthRsp(dotdeleted,:)=[];
        StimStrength=StrthRsp(:,1);     %Ψһ���ĸ��̼�ǿ��
        Respright=StrthRsp(:,4);        %�жϵ���ȷ��
        TrialCount=StrthRsp(:,3);       %ÿ���̼������ܵ��Դ���
        threshold(1,k)=fittingsigmoid(StimStrength,Respright,TrialCount); 
        figure_name=[subjectid '_' num2str(k)];
        saveas(gcf,[figure_name,'.jpg']);
        close all
        rawdata=[];
        StrthRsp=[];
    end
    save([subjectid ' Pure-tone Threshold'],"threshold")
end



%%�ú�������s���������
%%StimStrength ��ֵ
%%Respright     ��ȷ�жϵĸ���
%%TrialCount    ÿ���̼����trial��

function Threshold=fittingsigmoid(StimStrength,Respright,TrialCount)
    %function output=fittingsigmoid(StimStrength,Respright,TrialCount,subjectid,condition)
    
    
    %%%%%%%%%%%%%%��ϸ�˹�ֲ�����%%%%%%%%%%%%%%%%
    x = log10(min(StimStrength)):0.001:log10(max(StimStrength));
    %Threshold0= (log10(min(StrthRsp(:,1)))+ log10(max(StrthRsp(:,1))))/2;
    tempindex=find(TrialCount==max(TrialCount));
    Threshold0=log10(StimStrength(tempindex(1)));%��trial������ǿ����Ϊ�����ʼ����
    %Threshold0=log10(4);%��trial������ǿ����Ϊ�����ʼ����
    slope0=0.3;
    %amp0=3;
    beta0 = [Threshold0 slope0];
    fun = @(beta,x)(1+erf(((x-beta(1))/beta(2))/sqrt(2)))/2;                %%%chance level=0    ���
    %fun = @(beta,x)((1+erf(((x-beta(1))/beta(2))/sqrt(2)))/2)*0.5+0.5;     %%%chance level=0.5  ���
    %fun = @(beta,x)((1+erf(((x-beta(1))/beta(2))/sqrt(2)))/2)*0.75+0.25;   %%%chance level=0.25 ���
    
    %2014�汾matlab�������
    wnlm = fitnlm(log10(StimStrength),Respright,fun,beta0,'weight',TrialCount);
    r_square=wnlm.Rsquared.Ordinary;
    yfit = predict(wnlm,x');
    % SSR=sum(yfit-mean(Respright))*sum(yfit-mean(Respright));
    % SST=sum(Respright-mean(Respright))*sum(Respright-mean(Respright));
    % R2=SSR/SST;
    %%2011��֮ǰ�汾matlab�������
    % [wnlm r2] = nlinfit(log10(StimStrength),Respright,fun,beta0);
    % StimStrength
    % RespAccuracy
    % TrialCount
    % yfit=fun(wnlm,x');
    
    %%%%%%%%%%%%%%�� ͼ%%%%%%%%%%%%%%%%%%%%%%%%%
    %ͼ�ε���ɫ
    
    linecolorr=[0 0 1];
    dotcolorr=[0.5 0.5 0.8];
    
    
    %�����㡢��
    hold on
    plot(10.^x,yfit,'LineWidth',2,'Color',linecolorr)
    scatter(StimStrength,Respright,TrialCount*5,'ro','MarkerFaceColor',dotcolorr,'MarkerEdgeColor','none')
    ThreP=0.5;%1u1d
    ThreIndex=abs(yfit-ThreP)==min(abs(yfit-ThreP));
    Threshold=10^x(ThreIndex);%�ҵ���Ӧ����ֵ
    line([min(StimStrength) Threshold],[ThreP ThreP],'color',[0.5 0.5 0.5])
    line([Threshold Threshold],[ThreP 0],'color',dotcolorr)
    text(Threshold,0.05,num2str(Threshold,'%4.3g'));
    xlabel('Power_{Pure tone}/Power_{Noise} (dB)')
    ylabel('Probability')

end
