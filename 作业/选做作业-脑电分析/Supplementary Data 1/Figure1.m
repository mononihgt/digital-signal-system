% **************************************************
% Figure 1
% **************************************************

clear;clc

addpath ./tools

load Data/yeE1	%% EOG data in Experiment 1
x2=squeeze(mean(mean(ye(:,1,:,:),2),3));
v2=squeeze(mean(mean(ye(:,2,:,:),2),3));

load Data/脑电数据	%% EEG data in Experiment 1
yb = x;
yy=squeeze(mean(yb(:,:,:,:),3));
y2=squeeze(mean(mean(yb(:,48,:,:),2),3));clear y % channel Cz

fs=128/8;
t=0:1/fs:20;t=t-2;

%%
%EEG spectrum in experiment 1
Bv=1000;
wind=[t>=1 & t<11];
TD=sum(wind);
fy2=fft(y2(wind,:));  fy2=[2*abs(fy2)/TD].^2;
fx2=fft(x2(wind,:));  fx2=[2*abs(fx2)/TD].^2;
fv2=fft(v2(wind,:));  fv2=[2*abs(fv2)/TD].^2;
ff2=fft(yy(wind,:,:));ff2=[2*abs(ff2)/TD].^2;
[O,Afx2]=bootstrap_for_vector(Bv,0.05,@geomean,fx2');
[O,Afy2]=bootstrap_for_vector(Bv,0.05,@geomean,fy2');
[O,Afv2]=bootstrap_for_vector(Bv,0.05,@geomean,fv2');
f=1:size(fy2,1);f=f-1;f=f/size(fy2,1);f=f*fs;

figure;
subplot 111;hold on
plot_shade(f,geomean(fy2'),std(Afy2),[1 .5 .5])
plot(f,geomean(fy2,2),'r','linewidth',2)
xlim([.75 4.25]);ylim([0 .8]);
xlabel('Frequency (Hz)');ylabel('Power (μV^2)');

% subplot 312;hold on
% plot_shade(f,geomean(fx2'),std(Afx2),[.5 .5 1])
% plot(f,geomean(fx2,2),'b','linewidth',2)
% xlim([.75 4.25]);ylim([0 16])
% 
% subplot 313;hold on
% plot_shade(f,geomean(fv2'),std(Afv2),[.5 1 .5])
% plot(f,geomean(fv2,2),'color',[0 .7 0],'linewidth',2)
% xlim([.75 4.25]);ylim([0 21])


%%
%The significance of the spectral peaks in experiment 1
DIMs=[11 21 41];
Bv=1e4;
for dm=1:3
  DIM=DIMs(dm);
  dx2=[fx2(DIM,:); geomean(fx2([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:))];
  dy2=[fy2(DIM,:); geomean(fy2([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:))];
  dv2=[fv2(DIM,:); geomean(fv2([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:))];
  df2=[ff2(DIM,:,:); geomean(ff2([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,:))];
  if dm==1
    dxA=dx2;dyA=dy2;dvA=dv2;
  end
  
  [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dx2'); px2=[sum(A<=0)+1]/[Bv+1]; mx2=O;sx2=std(A);
  [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dy2'); py2=[sum(A<=0)+1]/[Bv+1]; my2=O;sy2=std(A);
  [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dv2'); pv2=[sum(A<=0)+1]/[Bv+1]; mv2=O;sv2=std(A);

  p(:,dm)=[py2 px2 pv2];
  m(:,dm)=[my2 mx2 mv2];
  s(:,dm)=[sy2 sx2 sv2];
end
[h, crit_p, p(1,:)]=fdr_bh(p(1,:));
[h, crit_p, p(2,:)]=fdr_bh(p(2,:));
[h, crit_p, p(3,:)]=fdr_bh(p(3,:));
for dm=1:3
  for id1=1:size(m,1)
    txtv{id1,dm}=[num2str(p(id1,dm)) ' (' num2str(m(id1,dm)) '/' num2str(s(id1,dm)) ')'];
  end
end
txtv

%%
%Correlation across the EEG, vertical EOG, and horizontal EOG
qx2=dxA(1,:)./dxA(2,:);
qy2=dyA(1,:)./dyA(2,:);
qv2=dvA(1,:)./dvA(2,:);
pv1=finv(1-0.05,2,8);
pv2=finv(1-0.005,2,8);

figure;
subplot 131;hold on
plot(qx2,qy2,'ro','linewidth',2)
[c,p]=corr(qx2',qy2');xlabel('vEOG');ylabel('EEG');title([c p])
set(gca,'xscale','log');set(gca,'yscale','log');axis([.01 200 .01 200]);set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);
[O,A]=bootstrap_for_vector(1e4,0.05,@c_corr,[qx2'+qy2'*1i]);[sum(A>0)+1]/10001

subplot 132;hold on
plot(qv2,qy2,'ro','linewidth',2)
[c,p]=corr(qv2',qy2');xlabel('hEOG');ylabel('EEG');title([c p])
set(gca,'xscale','log');set(gca,'yscale','log');axis([.01 200 .01 200]);set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);
[O,A]=bootstrap_for_vector(1e4,0.05,@c_corr,[qv2'+qy2'*1i]);[sum(A>0)+1]/10001

subplot 133;hold on
plot(qv2,qx2,'ro','linewidth',2)
[c,p]=corr(qv2',qx2');xlabel('hEOG');ylabel('vEOG');title([c p])
set(gca,'xscale','log');set(gca,'yscale','log');axis([.01 200 .01 200]);set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);
[O,A]=bootstrap_for_vector(1e4,0.05,@c_corr,[qv2'+qx2'*1i]);[sum(A>0)+1]/10001


