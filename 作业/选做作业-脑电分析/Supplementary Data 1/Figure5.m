% **************************************************
% Figure 5
% **************************************************

clear;clc

addpath .\tools

fs=16;
load Data\etkE4b  %% eyetracking data in Experiment 4b
y=squeeze(nanmean(ys(:,1,:,:,:,:),3));
y1(2*fs+[1:size(y,1)],:,:)=y;clear y
y=squeeze(nanmean(ys(:,2,:,:,:,:),3));
x1(2*fs+[1:size(y,1)],:,:)=y;clear y
y=squeeze(nanmean(ys(:,3,:,:,:,:),3));
z1(2*fs+[1:size(y,1)],:,:)=y;clear y
z1(isnan(z1))=nanmean(z1(:));

t=0:1/16:60;t=t(1:size(y1,1))-2;

%%
%blink, saccade, and pupil spectrum in 2 experimental conditions when eye is open in experiment 4b
tnow=[t>=7 & t<21];TD=sum(tnow);
fy0=fft(y1(tnow,:,:));fy0=[2*abs(fy0)/TD].^2;fy0=fy0+sqrt(mean(fy0(:).^2))/1000;
fx0=fft(x1(tnow,:,:));fx0=[2*abs(fx0)/TD].^2;fx0=fx0+sqrt(mean(fx0(:).^2))/1000;
fz0=fft(z1(tnow,:,:));fz0=[2*abs(fz0)/TD].^2;
f=1:size(fy0,1);f=f-1;f=16*f/size(fy0,1);

clrs(:,:,1)=[1 .7 .7];clrs(:,:,2)=[.7 .7 .7];
Bv=100;lnv={'r','k'};
figure;
for idx=1:2
  [O,Afx0]=bootstrap_for_vector(Bv,0.05,@geomean,fx0(:,:,idx)');
  [O,Afy0]=bootstrap_for_vector(Bv,0.05,@geomean,fy0(:,:,idx)');
  [O,Afz0]=bootstrap_for_vector(Bv,0.05,@geomean,fz0(:,:,idx)');
  
  subplot 311;hold on
  plot_shade(f,geomean(fy0(:,:,idx)'),std(Afy0),clrs(:,:,idx))
  plot(f,geomean(fy0(:,:,idx)'),lnv{idx},'linewidth',2)
  xlim([.25 2.25]);ylim([0 .005]);
  
  subplot 312;hold on
  plot_shade(f,geomean(fx0(:,:,idx)'),std(Afx0),clrs(:,:,idx))
  plot(f,geomean(fx0(:,:,idx)'),lnv{idx},'linewidth',2)
  xlim([.25 2.25]);ylim([0 2.2e-5]);
  
  subplot 313;hold on
  plot_shade(f,geomean(fz0(:,:,idx)'),std(Afz0),clrs(:,:,idx))
  plot(f,geomean(fz0(:,:,idx)'),lnv{idx},'linewidth',2)
  xlim([.25 2.25]);ylim([0 500]);
  set(gca,'ytick',0:200:1000)
end

%%
%The significance of blink, saccade, and pupil spectral peaks in 2 experimental conditions when eye is open in experiment 4b
DIMs=[8 15 29];
ffz=2;
for dimv=1:3
  for cc=1:2
    DIM=DIMs(dimv);
    dx0=[fx0(DIM,:,cc); geomean(fx0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];
    dy0=[fy0(DIM,:,cc); geomean(fy0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];
    dz0=[fz0(DIM,:,cc); geomean(fz0([[DIM-2]:[DIM-1],[DIM+1]:[DIM+0]],:,cc))];

    Bv=1e4;
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dx0'); px0=[sum(A<=0)+1]/[Bv+1];mx0=O;sx0=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dy0'); py0=[sum(A<=0)+1]/[Bv+1];my0=O;sy0=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dz0'); pz0=[sum(A<=0)+1]/[Bv+1];mz0=O;sz0=std(A);
  
    ps(:,cc,dimv)=[py0 px0 pz0];
    ms(:,cc,dimv)=[my0 mx0 mz0];
    ss(:,cc,dimv)=[sy0 sx0 sz0];
  end
end
for cc=1:2
  [h, crit_p, ps(1,cc,:)]=fdr_bh(squeeze(ps(1,cc,:)));
  [h, crit_p, ps(2,cc,:)]=fdr_bh(squeeze(ps(2,cc,:)));
  [h, crit_p, ps(3,cc,:)]=fdr_bh(squeeze(ps(3,cc,:)));
end

for cc=1:2
  for dm=1:3
    for id1=1:size(ms,1)
      txtv{id1,dm}=[num2str(ps(id1,cc,dm)) ' (' num2str(ms(id1,cc,dm)) '/' num2str(ss(id1,cc,dm)) ')'];
    end
  end
  txtv
end


%%
%The blink, saccade, and pupil response phase extracted by the Fourier analysisd in 2 experimental conditions when eye is open in experiment 4b
DIM=8;
fy0=fft(y1(tnow,:,:));fyF(:,:)=squeeze(fy0(DIM,:,:));
fx0=fft(x1(tnow,:,:));fxF(:,:)=squeeze(fx0(DIM,:,:));
fz0=fft(z1(tnow,:,:));fzF(:,:)=squeeze(fz0(DIM,:,:));
f=1:size(fy0,1);f=f-1;f=16*f/size(fy0,1);

figure;
subplot 311;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fyF(:,1)));
prcRGplot(angle(AA),95,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fyF(:,2)));
prcRGplot(angle(AA),95,'k');

subplot 312;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fxF(:,1)));
prcRGplot(angle(AA),95,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fxF(:,2)));
prcRGplot(angle(AA),95,'k');

subplot 313;hold on
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fzF(:,1)));
prcRGplot(angle(AA),95,'r');
[O1,AA] = bootstrap_for_vector(1e3,0.05,@mean,(fzF(:,2)));
prcRGplot(angle(AA),95,'k');

%%
%The blink, saccade, and pupil phase difference between conditions when eye is open in experiment 4b
figure;
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fxF);
px=prcRGplot(AA,95,'');
px2=prcRGplot(AA,99,'');
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fyF);
py=prcRGplot(AA,95,'');
py2=prcRGplot(AA,99,'');
[O1,AA] = bootstrap_for_vector(1e4,0.05,@meandf,fzF);
pz=prcRGplot(AA,95,'');
pz2=prcRGplot(AA,99,'');

A11=[py*180/pi;px*180/pi;pz*180/pi];
B11=[py2*180/pi;px2*180/pi;pz2*180/pi];
[A11(:,2) A11(:,[1 3]) B11(:,[1 3]) ]

%%
%Response waveform(0.25 ~ 0.75 Hz) for blink, saccade, and pupil when eye is open in experiment 4b
tg=[t>0 & t<21];
h=fir1(32,[.25 0.75]/8);
ly1(:,:,1)=filterBL(h,y1(:,:,1));lx1(:,:,1)=filterBL(h,x1(:,:,1));lz1(:,:,1)=filterBL(h,z1(:,:,1));
ly1(:,:,2)=filterBL(h,y1(:,:,2));lx1(:,:,2)=filterBL(h,x1(:,:,2));lz1(:,:,2)=filterBL(h,z1(:,:,2));
t0=t(1:size(ly1,1));

figure;
subplot 311;hold on
plot_shade(t0(tg),mean(ly1(tg,:,1)'),se(ly1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(ly1(tg,:,2)'),se(ly1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(ly1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(ly1(tg,:,2),2)),'k')
xlim([1 21]);ylim([-.1 1]*.3);set(gca,'xtick',0:2:22);grid on

subplot 312;hold on
plot_shade(t0(tg),mean(lx1(tg,:,1)'),se(lx1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(lx1(tg,:,2)'),se(lx1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lx1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(lx1(tg,:,2),2)),'k')
xlim([1 21]);ylim([0 .04]);set(gca,'xtick',0:2:22);grid on

subplot 313;hold on
plot_shade(t0(tg),mean(lz1(tg,:,1)'),se(lz1(tg,:,1)'),[1 .7 .7])
plot_shade(t0(tg),mean(lz1(tg,:,2)'),se(lz1(tg,:,2)'),[.7 .7 .7])
plot(t0(tg),squeeze(mean(lz1(tg,:,1),2)),'r')
plot(t0(tg),squeeze(mean(lz1(tg,:,2),2)),'k')
xlim([1 21]);ylim([-500 300]+2500);set(gca,'xtick',0:2:22);grid on

