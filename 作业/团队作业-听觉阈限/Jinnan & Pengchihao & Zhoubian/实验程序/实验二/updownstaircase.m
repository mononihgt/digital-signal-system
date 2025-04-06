%%Yongchun Cai ZJU
%%staircase procedure测阈值
%%staircasetype,阶梯的类型，字符串输入，1u1d、1u2d、1u3d
%%strengthmatrix,刺激强度矩阵
%%responsematrix，反应矩阵
%%minstrength，最小刺激
%%maxstrength，最长刺激 
%%根据最小刺激和最大刺激强度 确定步长大小，最小步长是大小刺激差别的1/12(单位log，即程序自动以log方式调节刺激强度)
%%initialstrength，刺激强度初始值

function stimulusstrength=updownstaircase(staircasetype,strengthmatrix,responsematrix,minstrength,maxstrength,initialstrength)
% 
% staircasetype='1u2d';
% strengthmatrix=tempstrength;
% responsematrix=tempresponse;
% minstrength=0.01;
% maxstrength=0.08;
% initialstrength=0.08;
minstep=abs(log10(maxstrength)-log10(minstrength))/12;          %将上下限分为12等分，也就是拟合阈值时会有最多13个数据点
max_div_min=8;
maxstep=minstep*max_div_min;                                    %最大步长为最小步长的8倍
steplist=minstep*[1 2 4 8 16];                                  %所有的step只能是minstep的2^n倍，2022-10-18
currentpoint=length(find(strengthmatrix~=0));                   %当前的trial序号（即最后一个强度不为0的trial）
%%
if currentpoint==0
   stimulusstrength=initialstrength; %stair的第一个trial
else
   switch staircasetype
      
      %%%%%%%%%%%%%%%1上1下阶梯法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      case '1u1d' 
          %%%%%%
          if responsematrix(currentpoint)==0
             stepsign=1; %按键错误，刺激强度上升
          elseif responsematrix(currentpoint)==1
             stepsign=-1; %按键正确，刺激强度下降
          end
          %%
          if currentpoint==1
             stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*maxstep);
          else
              %%%%%%%%%%%%本段代码用于寻找当前阶梯先前反转了几次，以确定当前的步长大小%%%%%%%%%%%%%%%
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
                     currentstep=abs(log10(strengthmatrix(currentpoint))-log10(strengthmatrix(currentpoint-1)));%当前的step
                        %%这里这样计算当前step的算法无法处理当刺激强度抵达上限/下限的时候，
                        %%因为此时当前刺激强度（上/下限）与上一个数据点的距离之间的距离不一定是之前的步长，而是小于等于先前的步长
        
                     %%%%%%%%%%2022-10-18%%%%%%%%%%%%
                     %%下面这段程序使得currentstep只能是minstep的1 2 4 8 16的倍数
                     [~,stepindex]=min(abs(steplist-currentstep));%找到currentstep与1 2 4 8 16被minstep最接近的step序号
                     currentstep=steplist(stepindex);
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %}

             if responsematrix(currentpoint-1)==responsematrix(currentpoint)
                stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
               elseif responsematrix(currentpoint-1)~=responsematrix(currentpoint)
                stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*max([currentstep/2 minstep]));
             end
          end
          
      %%%%%%%%%%%%%%%1上2下阶梯法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case '1u2d'
         if currentpoint==1
             if responsematrix(currentpoint)==0
                stepsign=1; %按键错误，刺激强度上升
             else
                stepsign=0; %按键正确，刺激强度不变
             end
             stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*maxstep);
         else
             %确定step的方向 
             if responsematrix(currentpoint)==0         %%一上
                  stepsign=1;
               elseif responsematrix(currentpoint-1)==1 && responsematrix(currentpoint)==1 && abs(log10(strengthmatrix(currentpoint-1))-log10(strengthmatrix(currentpoint)))<0.1*minstep %两次强度接近或相等时
                   stepsign=-1;                         %%二下
               else
                  stepsign=0;
             end
             
             %找最近的step大小
             laststep=0;
             kk=1;
             while abs(laststep)<0.1*minstep && (currentpoint-kk)>0 
               laststep=log10(strengthmatrix(currentpoint))-log10(strengthmatrix(currentpoint-kk));
               kk=kk+1;
             end
             
             %计算stimulusstrength
             if laststep==0 %stair刚开始，还没有step
                 currentstep=maxstep;
                 stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
             else
                %{
                             currentstep=abs(laststep);
                             %%%%%%%%%%2022-10-18%%%%%%%%%%%%
                             %%下面这段程序使得currentstep只能是minstep的1 2 4 8 16的倍数
                             [~,stepindex]=min(abs(steplist-currentstep));%找到currentstep与1 2 4 8 16被minstep最接近的step序号
                             currentstep=steplist(stepindex);
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %}
                  %%%%%%%%%%%%本段代码用于寻找当前阶梯先前反转了几次，以确定当前的步长大小%%%%%%%%%%%%%%%
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

                 if sign(laststep)==sign(stepsign) %与之前的step趋势一致，step不变
                    stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
                  elseif  sign(laststep)~=sign(stepsign) %与之前的step趋势相反，step减半，直到minstep
                    stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*max([currentstep/2 minstep]));
                 end       
             end
         end 
       
       %%%%%%%%%%%%%%%1上3下阶梯法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       case '1u3d' 
           last1point=currentpoint-1;
           last2point=currentpoint-2;
           if last1point<1||last2point<1            %%头3个试次内
               if responsematrix(currentpoint)==0
                 stepsign=1; %按键错误，刺激强度上升
                else
                 stepsign=0; %按键正确，刺激强度不变
               end
               stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*maxstep);
           else
               %确定step方向
               if responsematrix(currentpoint)==0
                   stepsign=1;
                elseif responsematrix(last1point)==1&&responsematrix(last2point)==1&&responsematrix(currentpoint)==1 && abs(log10(strengthmatrix(last2point))-log10(strengthmatrix(currentpoint)))<0.1*minstep %三次强度接近或相等时
                   stepsign=-1;
                else
                   stepsign=0;
               end
         
               %找最近的step大小
               laststep=0;
               kk=1;
               while abs(laststep)<0.1*minstep && (currentpoint-kk)>0
                 laststep=log10(strengthmatrix(currentpoint))-log10(strengthmatrix(currentpoint-kk));
                 kk=kk+1;
               end
               
               %计算stimulusstrength
               if laststep==0 %stair刚开始，还没有step
                  currentstep=maxstep;
                  stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
                else
                  %{
                    currentstep=abs(laststep);
                  %%%%%%%%%%2022-10-18%%%%%%%%%%%%
                  %%下面这段程序使得currentstep只能是minstep的1 2 4 8 16的倍数
                  [~,stepindex]=min(abs(steplist-currentstep));%找到currentstep与1 2 4 8 16被minstep最接近的step序号
                  currentstep=steplist(stepindex);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %}
                  %%%%%%%%%%%%本段代码用于寻找当前阶梯先前反转了几次，以确定当前的步长大小%%%%%%%%%%%%%%%
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
                  if sign(laststep)==sign(stepsign) %与之前的step趋势一致，step不变
                     stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*currentstep);
                   elseif  sign(laststep)~=sign(stepsign) %与之前的step趋势相反，step减半，直到minstep
                     stimulusstrength=10^(log10(strengthmatrix(currentpoint))+stepsign*max([currentstep/2 minstep]));
                  end       
               end
               
           end
       %%%%%%%%%%%%%%%%%%%%%%%
   end
          
end


%%%限制刺激范围在minstrength到maxstrength之内
stimulusstrength=max([minstrength stimulusstrength]);
stimulusstrength=min([maxstrength stimulusstrength]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    