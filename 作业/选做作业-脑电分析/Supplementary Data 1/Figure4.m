% **************************************************
% Figure 4
% **************************************************

clear;clc

addpath .\tools

fs=16;
load Data\yeE4	%% EEG data in Experiment 4
ye(abs(ye)>1e3)=0;
x1=squeeze(mean(mean(ye(:,1:2,:,:,:),2),3));
z1=squeeze(mean(mean(ye(:,3,:,:,:),2),3));

load Data\etkE4	%% eyetracking data in Experiment 4
y1=squeeze(nanmean(ys(:,1,:,:,:,:),3));
u1=squeeze(nanmean(ys(:,2,:,:,:,:),3));
w1=squeeze(nanmean(ys(:,3,:,:,:,:),3));
ffz=1;

y1(2*fs+[1:size(ys,1)],:,:)=y1;y1(size(x1,1),:,:)=0;y1(:,:,4)=0;
u1(2*fs+[1:size(ys,1)],:,:)=u1;u1(size(x1,1),:,:)=0;u1(:,:,4)=0;
w1(2*fs+[1:size(ys,1)],:,:)=w1;w1(size(x1,1),:,:)=0;w1(w1==0)=median(w1(:));w1(:,:,4)=0;

t=0:1/16:20;t=t(1:size(y1,1))-2;

%%

tnow=[t>=4 & t<12];TD=sum(tnow);
fy0=fft(y1(tnow,:,:));fy0=[2*abs(fy0)/TD].^2;fy0=fy0+sqrt(mean(fy0(:).^2))/1000;
fx0=fft(x1(tnow,:,:));fx0=[2*abs(fx0)/TD].^2;
fz0=fft(z1(tnow,:,:));fz0=[2*abs(fz0)/TD].^2;
fu0=fft(u1(tnow,:,:));fu0=[2*abs(fu0)/TD].^2;fu0=fu0+sqrt(mean(fu0(:).^2))/1000;
fw0=fft(w1(tnow,:,:));fw0=[2*abs(fw0)/TD].^2;
f=1:size(fy0,1);f=f-1;f=16*f/size(fy0,1);

Bv=100;lnv={'r','k','r','k'};clrs(:,:,1)=[1 .7 .7];clrs(:,:,2)=[.7 .7 .7];
figure;
%vertical EOG, horizontal EOG, blink, saccade, and pupil spectrum in 2 experimental conditions when eye is open in experiment 4
for idx=1:2
  [O,Afx0]=bootstrap_for_vector(Bv,0.05,@geomean,fx0(:,:,idx)');
  [O,Afy0]=bootstrap_for_vector(Bv,0.05,@geomean,fy0(:,:,idx)');
  [O,Afz0]=bootstrap_for_vector(Bv,0.05,@geomean,fz0(:,:,idx)');
  [O,Afu0]=bootstrap_for_vector(Bv,0.05,@geomean,fu0(:,:,idx)');
  [O,Afw0]=bootstrap_for_vector(Bv,0.05,@geomean,fw0(:,:,idx)');
  
  subplot 711;hold on;plot_shade(f,geomean(fx0(:,:,idx)'),std(Afx0),clrs(:,:,idx))
  plot(f,geomean(fx0(:,:,idx)'),lnv{idx},'linewidth',2);xlim([.75 4.25]);ylim([0 75]);
  set(gca,'ytick',0:20:100);
  
  subplot 712;hold on;plot_shade(f,geomean(fz0(:,:,idx)'),std(Afz0),clrs(:,:,idx))
  plot(f,geomean(fz0(:,:,idx)'),lnv{idx},'linewidth',2);xlim([.75 4.25]);ylim([0 2.5]);
  set(gca,'ytick',0:3);
  
  subplot 713;hold on;plot_shade(f,geomean(fy0(:,:,idx)'),std(Afy0),clrs(:,:,idx))
  plot(f,geomean(fy0(:,:,idx)'),lnv{idx},'linewidth',2);xlim([.75 4.25]);%ylim([0 .0005]);
  set(gca,'ytick',[0:2:6]*1e-4);
  
  subplot 714;hold on;plot_shade(f,geomean(fu0(:,:,idx)'),std(Afu0),clrs(:,:,idx))
  plot(f,geomean(fu0(:,:,idx)'),lnv{idx},'linewidth',2);xlim([.75 4.25]);ylim([0 .00002]);
  
  subplot 715;hold on;plot_shade(f([f>.2 & f<10]),geomean(fw0([f>.2 & f<10],:,idx)'),std(Afw0(:,[f>.2 & f<10])),clrs(:,:,idx))
  plot(f([f>.2 & f<10]),geomean(fw0([f>.2 & f<10],:,idx)'),lnv{idx},'linewidth',2);xlim([.75 4.25]);ylim([0 70]);
  set(gca,'ytick',0:20:60);
end
%vertical EOG, horizontal EOG spectrum in 2 experimental conditions when eye is close in experiment 4
for idx=3:4
  [O,Afx0]=bootstrap_for_vector(Bv,0.05,@geomean,fx0(:,:,idx)');
  [O,Afz0]=bootstrap_for_vector(Bv,0.05,@geomean,fz0(:,:,idx)');
  
  subplot 716;hold on;plot_shade(f,geomean(fx0(:,:,idx)'),std(Afx0),clrs(:,:,idx-2))
  plot(f,geomean(fx0(:,:,idx)'),lnv{idx},'linewidth',2);xlim([.75 4.25]);ylim([0 25]);
  set(gca,'ytick',0:10:100);
  
  subplot 717;hold on;plot_shade(f,geomean(fz0(:,:,idx)'),std(Afz0),clrs(:,:,idx-2))
  plot(f,geomean(fz0(:,:,idx)'),lnv{idx},'linewidth',2);xlim([.75 4.25]);ylim([0 7]);
  set(gca,'ytick',0:2:7);
end

%% peak significance
%The significance of vertical EOG, horizontal EOG, blink, saccade, and pupil spectral peaks in 2 experimental conditions when eye is open in experiment 4
%The significance of vertical EOG, horizontal EOG spectral peaks in 2 experimental conditions when eye is close in experiment 4
DIMs=[9 17 33];
ffz=2;
for dimv=1:3
  for cc=1:4
    DIM=DIMs(dimv);
    dx0=[fx0(DIM,:,cc); geomean(fx0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];
    dy0=[fy0(DIM,:,cc); geomean(fy0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];
    dz0=[fz0(DIM,:,cc); geomean(fz0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];
    du0=[fu0(DIM,:,cc); geomean(fu0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];
    dw0=[fw0(DIM,:,cc); geomean(fw0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];

    Bv=1e4;
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dx0'); px0=[sum(A<=0)+1]/[Bv+1];mx0=O;sx0=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dy0'); py0=[sum(A<=0)+1]/[Bv+1];my0=O;sy0=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dz0'); pz0=[sum(A<=0)+1]/[Bv+1];mz0=O;sz0=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,du0'); pu0=[sum(A<=0)+1]/[Bv+1];mu0=O;su0=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dw0'); pw0=[sum(A<=0)+1]/[Bv+1];mw0=O;sw0=std(A);
  
    ps(:,cc,dimv)=[px0 pz0 py0 pu0 pw0];
    ms(:,cc,dimv)=[mx0 mz0 my0 mu0 mw0];
    ss(:,cc,dimv)=[sx0 sz0 sy0 su0 sw0];
  end
end
for cc=1:4
  [h, crit_p, ps(1,cc,:)]=fdr_bh(squeeze(ps(1,cc,:)));
  [h, crit_p, ps(2,cc,:)]=fdr_bh(squeeze(ps(2,cc,:)));
  [h, crit_p, ps(3,cc,:)]=fdr_bh(squeeze(ps(3,cc,:)));
end

for cc=1:4
  for dm=1:3
    for id1=1:size(ms,1)
      txtv{id1,dm}=[num2str(ps(id1,cc,dm)) ' (' num2str(ms(id1,cc,dm)) '/' num2str(ss(id1,cc,dm)) ')'];
    end
  end
  txtv
end


%%
DIM=9;
tnow=[t>=4 & t<12];
fy0=fft(y1(tnow,:,:));fyF(:,:)=squeeze(fy0(DIM,:,:));
fx0=fft(x1(tnow,:,:));fxF(:,:)=squeeze(fx0(DIM,:,:));
fz0=fft(z1(tnow,:,:));fzF(:,:)=squeeze(fz0(DIM,:,:));
fu0=fft(u1(tnow,:,:));fuF(:,:)=squeeze(fu0(DIM,:,:));
fw0=fft(w1(tnow,:,:));fwF(:,:)=squeeze(fw0(DIM,:,:));
f=1:size(fy0,1);f=f-1;f=16*f/size(fy0,1);

figure;
%The vertical EOG, horizontal EOG, blink, saccade, and pupil response phase extracted by the Fourier analysisd in 2 experimental conditions when eye is open in experiment 4
subplot 711;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fxF(:,1)));prcRGplot(angle(AA),90,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fxF(:,2)));prcRGplot(angle(AA),90,'k');
subplot 712;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fzF(:,1)));prcRGplot(angle(AA),90,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fzF(:,2)));prcRGplot(angle(AA),90,'k');
subplot 713;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fyF(:,1)));prcRGplot(angle(AA),90,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fyF(:,2)));prcRGplot(angle(AA),90,'k');
subplot 714;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fuF(:,1)));prcRGplot(angle(AA),90,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fuF(:,2)));prcRGplot(angle(AA),90,'k');
subplot 715;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fwF(:,1)));prcRGplot(angle(AA),90,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fwF(:,2)));prcRGplot(angle(AA),90,'k');

%The vertical EOG, horizontal EOG response phase extracted by the Fourier analysisd in 2 experimental conditions when eye is close in experiment 4
subplot 716;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fxF(:,3)));prcRGplot(angle(AA),90,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fxF(:,4)));prcRGplot(angle(AA),90,'k');
subplot 717;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fzF(:,3)));prcRGplot(angle(AA),90,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fzF(:,4)));prcRGplot(angle(AA),90,'k');

%%
%The significance of vertical EOG, horizontal EOG, blink, saccade, and pupil phase difference between conditions when eye is open in experiment 4
figure;
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fxF(:,1:2));
px=prcRGplot(AA,95,'');
px2=prcRGplot(AA,99,'');
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fzF(:,1:2));
pz=prcRGplot(AA,95,'');
pz2=prcRGplot(AA,99,'');
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fyF(:,1:2));
py=prcRGplot(AA,95,'');
py2=prcRGplot(AA,99,'');
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fuF(:,1:2));
pu=prcRGplot(AA,95,'');
pu2=prcRGplot(AA,99,'');
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fwF(:,1:2));
pw=prcRGplot(AA,95,'');
pw2=prcRGplot(AA,99,'');

A11=[px*180/pi;pz*180/pi;py*180/pi;pu*180/pi;pw*180/pi];
B11=[px2*180/pi;pz2*180/pi;py2*180/pi;pu2*180/pi;pw2*180/pi];
[A11(:,2) A11(:,[1 3]) B11(:,[1 3]) ]

%%
%The significance of vertical EOG, horizontal EOG phase difference between conditions when eye is close in experiment 4
figure;
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fxF(:,3:4));
px=prcRGplot(AA,95,'');
px2=prcRGplot(AA,99,'');
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fzF(:,3:4));
pz=prcRGplot(AA,95,'');
pz2=prcRGplot(AA,99,'');

A11=[px*180/pi;pz*180/pi];
B11=[px2*180/pi;pz2*180/pi];
[A11(:,2) A11(:,[1 3]) B11(:,[1 3]) ]


%%
h=fir1(24,[.75 1.25]/8);
ly1(:,:,1)=filterBL(h,y1(:,:,1));lx1(:,:,1)=filterB(h,x1(:,:,1));lz1(:,:,1)=filterB(h,z1(:,:,1));
ly1(:,:,2)=filterBL(h,y1(:,:,2));lx1(:,:,2)=filterB(h,x1(:,:,2));lz1(:,:,2)=filterB(h,z1(:,:,2));
ly1(:,:,3)=filterBL(h,y1(:,:,3));lx1(:,:,3)=filterB(h,x1(:,:,3));lz1(:,:,3)=filterB(h,z1(:,:,3));
ly1(:,:,4)=filterBL(h,y1(:,:,4));lx1(:,:,4)=filterB(h,x1(:,:,4));lz1(:,:,4)=filterB(h,z1(:,:,4));
lu1(:,:,1)=filterBL(h,u1(:,:,1));lw1(:,:,1)=filterBL(h,w1(:,:,1));
lu1(:,:,2)=filterBL(h,u1(:,:,2));lw1(:,:,2)=filterBL(h,w1(:,:,2));

t0=t(1:size(ly1,1));
tg=[t0>0 & t0<12];
figure;
%Response waveform(0.75 ~ 1.25 Hz) for vertical EOG, horizontal EOG, blink, saccade, and pupil when eye is open in experiment 4
subplot 711;hold on
plot_shade(t0(tg),mean(lx1(tg,:,1)'),se(lx1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(lx1(tg,:,2)'),se(lx1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lx1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(lx1(tg,:,2),2)),'k')
xlim([0 12]);ylim([-25 25]);set(gca,'xtick',0:12);grid on

subplot 712;hold on
plot_shade(t0(tg),mean(lz1(tg,:,1)'),se(lz1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(lz1(tg,:,2)'),se(lz1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lz1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(lz1(tg,:,2),2)),'k')
xlim([0 12]);ylim([-5 5]);set(gca,'xtick',0:12);grid on

subplot 713;hold on
plot_shade(t0(tg),mean(ly1(tg,:,1)'),se(ly1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(ly1(tg,:,2)'),se(ly1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(ly1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(ly1(tg,:,2),2)),'k')
xlim([0 12]);ylim([0 1]*.2);set(gca,'xtick',0:12);grid on

subplot 714;hold on
plot_shade(t0(tg),mean(lu1(tg,:,1)'),se(lu1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(lu1(tg,:,2)'),se(lu1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lu1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(lu1(tg,:,2),2)),'k')
xlim([0 12]);ylim([0 40]*1e-3);set(gca,'xtick',0:12);grid on

subplot 715;hold on
plot_shade(t0(tg),mean(lw1(tg,:,1)'),se(lw1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(lw1(tg,:,2)'),se(lw1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lw1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(lw1(tg,:,2),2)),'k')
xlim([0 12]);ylim([1.8e3 2.3e3]);set(gca,'xtick',0:12);grid on

%Response waveform(0.75 ~ 1.25 Hz) for vertical EOG, horizontal EOG when eye is close in experiment 4
subplot 716;hold on
plot_shade(t0(tg),mean(lx1(tg,:,3)'),se(lx1(tg,:,3)'),[1 .7 .7])
plot_shade(t0(tg),mean(lx1(tg,:,4)'),se(lx1(tg,:,4)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lx1(tg,:,3),2)),'r')
plot(t0(tg),squeeze(mean(lx1(tg,:,4),2)),'k')
xlim([0 12]);ylim([-10 10]);set(gca,'xtick',0:12);grid on

subplot 717;hold on
plot_shade(t0(tg),mean(lz1(tg,:,3)'),se(lz1(tg,:,3)'),[1 .7 .7])
plot_shade(t0(tg),mean(lz1(tg,:,4)'),se(lz1(tg,:,4)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lz1(tg,:,3),2)),'r')
plot(t0(tg),squeeze(mean(lz1(tg,:,4),2)),'k')
xlim([0 12]);ylim([-5 5]);set(gca,'xtick',0:12);grid on
