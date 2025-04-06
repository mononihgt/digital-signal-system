%%Yongchun Cai ZJU
%%staircase procedure����ֵ
%%staircasetype,���ݵ����ͣ��ַ������룬1u1d��1u2d��1u3d
%%strengthmatrix,�̼�ǿ�Ⱦ���
%%responsematrix����Ӧ����
%%minstrength����С�̼�
%%maxstrength����̼� 
%%������С�̼������̼�ǿ�� ȷ��������С����С�����Ǵ�С�̼�����1/12(��λlog���������Զ���log��ʽ���ڴ̼�ǿ��)
%%initialstrength���̼�ǿ�ȳ�ʼֵ

function stimulusstrength=updownstaircase(staircasetype,strengthmatrix,responsematrix,minstrength,maxstrength,initialstrength)
% 
% staircasetype='1u2d';
% strengthmatrix=tempstrength;
% responsematrix=tempresponse;
% minstrength=0.01;
% maxstrength=0.08;
% initialstrength=0.08;
minstep=abs(log10(maxstrength)-log10(minstrength))/12;          %�������޷�Ϊ12�ȷ֣�Ҳ���������ֵʱ�������13�����ݵ�
max_div_min=8;
maxstep=minstep*max_div_min;                                    %��󲽳�Ϊ��С������8��
steplist=minstep*[1 2 4 8 16];                                  %���е�stepֻ����minstep��2^n����2022-10-18
currentpoint=length(find(strengthmatrix~=0));                   %��ǰ��trial��ţ������һ��ǿ�Ȳ�Ϊ0��trial��
%%
if currentpoint==0
   stimulusstrength=initialstrength; %stair�ĵ�һ��trial
else
   switch staircasetype
      
      %%%%%%%%%%%%%%%1��1�½��ݷ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      case '1u1d' 
          %%%%%%
          if responsematrix(currentpoint)==0
             stepsign=1; %�������󣬴̼�ǿ������
          elseif responsematrix(currentpoint)==1
             stepsign=-1; %������ȷ���̼�ǿ���½�
          end
          %%
          if currentpoint==1
             stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*maxstep);
          else
              %%%%%%%%%%%%���δ�������Ѱ�ҵ�ǰ������ǰ��ת�˼��Σ���ȷ����ǰ�Ĳ�����С%%%%%%%%%%%%%%%
              kk=2;
              last_diff=0;
              reverse_time=0;
              while kk<=currentpoint && reverse_time<log2(max_div_min)
%                   last_diff
                  current_diff=log10(strengthmatrix(kk))-log10(strengthmatrix(kk-1));
                  if current_diff~=0 && abs(current_diff)>0.1*minstep
                      
                      if sign(last_diff)*sign(current_diff)<0
                        
                          reverse_time=reverse_time+1;
                      end

                      last_diff=current_diff;
                  end
                  kk=kk+1;
              end
              currentstep=steplist(log2(max_div_min)+1-reverse_time);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %{
                     currentstep=abs(log10(strengthmatrix(currentpoint))-log10(strengthmatrix(currentpoint-1)));%��ǰ��step
                        %%�����������㵱ǰstep���㷨�޷������̼�ǿ�ȵִ�����/���޵�ʱ��
                        %%��Ϊ��ʱ��ǰ�̼�ǿ�ȣ���/���ޣ�����һ�����ݵ�ľ���֮��ľ��벻һ����֮ǰ�Ĳ���������С�ڵ�����ǰ�Ĳ���
        
                     %%%%%%%%%%2022-10-18%%%%%%%%%%%%
                     %%������γ���ʹ��currentstepֻ����minstep��1 2 4 8 16�ı���
                     [~,stepindex]=min(abs(steplist-currentstep));%�ҵ�currentstep��1 2 4 8 16��minstep��ӽ���step���
                     currentstep=steplist(stepindex);
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %}

             if responsematrix(currentpoint-1)==responsematrix(currentpoint)
                stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
               elseif responsematrix(currentpoint-1)~=responsematrix(currentpoint)
                stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*max([currentstep/2 minstep]));
             end
          end
          
      %%%%%%%%%%%%%%%1��2�½��ݷ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case '1u2d'
         if currentpoint==1
             if responsematrix(currentpoint)==0
                stepsign=1; %�������󣬴̼�ǿ������
             else
                stepsign=0; %������ȷ���̼�ǿ�Ȳ���
             end
             stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*maxstep);
         else
             %ȷ��step�ķ��� 
             if responsematrix(currentpoint)==0         %%һ��
                  stepsign=1;
               elseif responsematrix(currentpoint-1)==1 && responsematrix(currentpoint)==1 && abs(log10(strengthmatrix(currentpoint-1))-log10(strengthmatrix(currentpoint)))<0.1*minstep %����ǿ�Ƚӽ������ʱ
                   stepsign=-1;                         %%����
               else
                  stepsign=0;
             end
             
             %�������step��С
             laststep=0;
             kk=1;
             while abs(laststep)<0.1*minstep && (currentpoint-kk)>0 
               laststep=log10(strengthmatrix(currentpoint))-log10(strengthmatrix(currentpoint-kk));
               kk=kk+1;
             end
             
             %����stimulusstrength
             if laststep==0 %stair�տ�ʼ����û��step
                 currentstep=maxstep;
                 stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
             else
                %{
                             currentstep=abs(laststep);
                             %%%%%%%%%%2022-10-18%%%%%%%%%%%%
                             %%������γ���ʹ��currentstepֻ����minstep��1 2 4 8 16�ı���
                             [~,stepindex]=min(abs(steplist-currentstep));%�ҵ�currentstep��1 2 4 8 16��minstep��ӽ���step���
                             currentstep=steplist(stepindex);
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %}
                  %%%%%%%%%%%%���δ�������Ѱ�ҵ�ǰ������ǰ��ת�˼��Σ���ȷ����ǰ�Ĳ�����С%%%%%%%%%%%%%%%
                  kk=2;
                  last_diff=0;
                  reverse_time=0;
                  while kk<=currentpoint && reverse_time<log2(max_div_min)
                      current_diff=log10(strengthmatrix(kk))-log10(strengthmatrix(kk-1));
                      if current_diff~=0 && abs(current_diff)>0.1*minstep
                          if sign(last_diff)*sign(current_diff)<0
                              reverse_time=reverse_time+1;
                          end
                          last_diff=current_diff;
                      end
                      kk=kk+1;
                  end
                  currentstep=steplist(log2(max_div_min)+1-reverse_time);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

                 if sign(laststep)==sign(stepsign) %��֮ǰ��step����һ�£�step����
                    stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
                  elseif  sign(laststep)~=sign(stepsign) %��֮ǰ��step�����෴��step���룬ֱ��minstep
                    stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*max([currentstep/2 minstep]));
                 end       
             end
         end 
       
       %%%%%%%%%%%%%%%1��3�½��ݷ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       case '1u3d' 
           last1point=currentpoint-1;
           last2point=currentpoint-2;
           if last1point<1||last2point<1            %%ͷ3���Դ���
               if responsematrix(currentpoint)==0
                 stepsign=1; %�������󣬴̼�ǿ������
                else
                 stepsign=0; %������ȷ���̼�ǿ�Ȳ���
               end
               stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*maxstep);
           else
               %ȷ��step����
               if responsematrix(currentpoint)==0
                   stepsign=1;
                elseif responsematrix(last1point)==1&&responsematrix(last2point)==1&&responsematrix(currentpoint)==1 && abs(log10(strengthmatrix(last2point))-log10(strengthmatrix(currentpoint)))<0.1*minstep %����ǿ�Ƚӽ������ʱ
                   stepsign=-1;
                else
                   stepsign=0;
               end
         
               %�������step��С
               laststep=0;
               kk=1;
               while abs(laststep)<0.1*minstep && (currentpoint-kk)>0
                 laststep=log10(strengthmatrix(currentpoint))-log10(strengthmatrix(currentpoint-kk));
                 kk=kk+1;
               end
               
               %����stimulusstrength
               if laststep==0 %stair�տ�ʼ����û��step
                  currentstep=maxstep;
                  stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
                else
                  %{
                    currentstep=abs(laststep);
                  %%%%%%%%%%2022-10-18%%%%%%%%%%%%
                  %%������γ���ʹ��currentstepֻ����minstep��1 2 4 8 16�ı���
                  [~,stepindex]=min(abs(steplist-currentstep));%�ҵ�currentstep��1 2 4 8 16��minstep��ӽ���step���
                  currentstep=steplist(stepindex);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %}
                  %%%%%%%%%%%%���δ�������Ѱ�ҵ�ǰ������ǰ��ת�˼��Σ���ȷ����ǰ�Ĳ�����С%%%%%%%%%%%%%%%
                  kk=2;
                  last_diff=0;
                  reverse_time=0;
                  while kk<=currentpoint && reverse_time<log2(max_div_min)
                      current_diff=log10(strengthmatrix(kk))-log10(strengthmatrix(kk-1));
                      if current_diff~=0 && abs(current_diff)>0.1*minstep
                          if sign(last_diff)*sign(current_diff)<0
                              reverse_time=reverse_time+1;
                          end
                          last_diff=current_diff;
                      end
                      kk=kk+1;
                  end
                  currentstep=steplist(log2(max_div_min)+1-reverse_time);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                  if sign(laststep)==sign(stepsign) %��֮ǰ��step����һ�£�step����
                     stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
                   elseif  sign(laststep)~=sign(stepsign) %��֮ǰ��step�����෴��step���룬ֱ��minstep
                     stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*max([currentstep/2 minstep]));
                  end       
               end
               
           end
       %%%%%%%%%%%%%%%%%%%%%%%
   end
          
end


%%%���ƴ̼���Χ��minstrength��maxstrength֮��
stimulusstrength=max([minstrength stimulusstrength]);
stimulusstrength=min([maxstrength stimulusstrength]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    