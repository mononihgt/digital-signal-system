% **************************************************
% Figure 2
% **************************************************

clear;

addpath .\tools

fs=128/8;
load Data\yeE2	%% EOG data in Experiment 2
x4=squeeze(mean(mean(ye(:,1,:,:,:),2),3));
w4=squeeze(mean(mean(ye(:,2,:,:,:),2),3));

load Data\ybE2	%% EEG data in Experiment 2
y4=squeeze(mean(mean(yb(:,5,:,:,:),2),3));clear y % channel Cz

load Data\etkE2;	%% eyetracking data in Experiment 2
y=squeeze(nanmean(ys(:,1,:,:,:,:),3));
v4(2*fs+[1:size(y,1)],:,[1:2 4:5])=y;clear y
y=squeeze(nanmean(ys(:,2,:,:,:,:),3));
u4(2*fs+[1:size(y,1)],:,[1:2 4:5])=y;clear y
y=squeeze(nanmean(ys(:,3,:,:,:,:),3));
z4(2*fs+[1:size(y,1)],:,[1:2 4:5])=y;clear y
for id3=1:size(z4,3)
  for id2=1:size(z4,2)
    zz4=z4(:,id2,id3);
    zz4(isnan(zz4))=nanmean(zz4);
    z4(:,id2,id3)=zz4;
  end
end
ffz=1;

%%
%Response spectra in 4 experimental conditions (columns) and for 6 neural/ocular measures (rows) in experiment 2
t=0:1/fs:20;t=t(1:size(x4))-2;
cdts={3,4,[1:2],5};lst='-';

figure;
for cc=1:4
  cdt=cdts{cc};
  Bv=100;
  wind=[t>=1 & t<11];
  TD=sum(wind);
  fy4=fft(mean(y4(wind,:,cdt),3));fy4=[2*abs(fy4)/TD].^2;
  fx4=fft(mean(x4(wind,:,cdt),3));fx4=[2*abs(fx4)/TD].^2;
  fv4=fft(mean(v4(wind,:,cdt),3));fv4=[2*abs(fv4)/TD].^2;fv4=fv4+sqrt(mean(fv4(:).^2))/1000;
  fw4=fft(mean(w4(wind,:,cdt),3));fw4=[2*abs(fw4)/TD].^2;
  fu4=fft(mean(u4(wind,:,cdt),3));fu4=[2*abs(fu4)/TD].^2;fu4=fu4+sqrt(mean(fu4(:).^2))/1000;
  fz4=fft(mean(z4(wind,:,cdt),3));fz4=[2*abs(fz4)/TD].^2;
  [O,Afx4]=bootstrap_for_vector(Bv,0.05,@geomean,fx4');
  [O,Afy4]=bootstrap_for_vector(Bv,0.05,@geomean,fy4');
  try
    [O,Afv4]=bootstrap_for_vector(Bv,0.05,@geomean,fv4');
  catch Afv4=Afy4*0;end
  [O,Afw4]=bootstrap_for_vector(Bv,0.05,@geomean,fw4');
  [O,Afu4]=bootstrap_for_vector(Bv,0.05,@geomean,fu4');
  [O,Afz4]=bootstrap_for_vector(Bv,0.05,@geomean,fz4');
  
  f=1:size(fy4,1);f=f-1;f=f/size(fy4,1);f=f*fs;
  subplot(6,4,cc);hold on
  plot_shade(f(2:end),geomean(fy4(2:end,:)'),std(Afy4(:,2:end)),[1 .5 .5])
  plot(f,geomean(fy4,2),'r','linewidth',2,'linestyle',lst)
  xlim([.7 4.25]);ylim([0 .8])
  set(gca,'ytick',0:.2:.8);set(gca,'xtick',[1 2 4]);
  
  subplot(6,4,4+cc);hold on
  plot_shade(f(2:end),geomean(fx4(2:end,:)'),std(Afx4(:,2:end)),[.5 .5 1])
  plot(f,geomean(fx4,2),'b','linewidth',2,'linestyle',lst)
  xlim([.7 4.25]);ylim([0 64])
  set(gca,'ytick',0:20:120);set(gca,'xtick',[1 2 4]);
  
  subplot(6,4,8+cc);hold on
  plot_shade(f(2:end),geomean(fw4(2:end,:)'),std(Afw4(:,2:end)),[.5 1 .5])
  plot(f,geomean(fw4,2),'color',[0 .7 0],'linewidth',2,'linestyle',lst)
  xlim([.7 4.25]);ylim([0 3.5])
  set(gca,'ytick',0:1:20);set(gca,'xtick',[1 2 4]);
  
  subplot(6,4,12+cc);hold on;
  plot_shade(f(2:end),geomean(fv4(2:end,:)'),std(Afv4(:,2:end)),[1 .8 .3])
  plot(f,geomean(fv4,2),'color',[1 .6 0],'linewidth',2,'linestyle',lst)
  xlim([.7 4.25]);ylim([0 0.00085])
  set(gca,'ytick',[0:3:10]*1e-4);set(gca,'xtick',[1 2 4]);
  
  subplot(6,4,16+cc);hold on;
  plot_shade(f(2:end),geomean(fu4(2:end,:)'),std(Afu4(:,2:end)),[.8 .5 1])
  plot(f,geomean(fu4,2),'color',[.7 .2 1],'linewidth',2,'linestyle',lst)
  xlim([.7 4.25]);ylim([0 1.5e-5])
  set(gca,'ytick',[0:.5:10]*1e-5);set(gca,'xtick',[1 2 4]);
  
  subplot(6,4,20+cc);hold on;
  plot_shade(f(2:end),geomean(fz4(2:end,:)'),std(Afz4(:,2:end)),[.5 1 1])
  plot(f,geomean(fz4,2),'color',[0 .6 1],'linewidth',2,'linestyle',lst)
  xlim([.7 4.25]);ylim([0 400])
  set(gca,'ytick',[0:3:10]*50);set(gca,'xtick',[1 2 4]);
end
subplot(6,4,5);ylim([0 6]);set(gca,'ytick',0:2:20);
subplot(6,4,9);ylim([0 11]);set(gca,'ytick',0:3:20);
subplot(6,4,13);cla;axis off
subplot(6,4,17);cla;axis off
subplot(6,4,21);cla;axis off

%%
%The significance of the spectral peaks in experiment 2
cdts={3,4,[1:2],5};
DIMs=[11 21 41];
ffz=2;
for dimv=1:3
  for cc=1:4
    cdt=cdts{cc};
    fy4=fft(mean(y4(wind,:,cdt),3));fy4=[2*abs(fy4)/TD].^2;
    fx4=fft(mean(x4(wind,:,cdt),3));fx4=[2*abs(fx4)/TD].^2;
    fw4=fft(mean(w4(wind,:,cdt),3));fw4=[2*abs(fw4)/TD].^2;
    fv4=fft(mean(v4(wind,:,cdt),3));fv4=[2*abs(fv4)/TD].^2;fv4=fv4+sqrt(mean(fv4(:).^2))/1000;
    fu4=fft(mean(u4(wind,:,cdt),3));fu4=[2*abs(fu4)/TD].^2;fu4=fu4+sqrt(mean(fu4(:).^2))/1000;
    fz4=fft(mean(z4(wind,:,cdt),3));fz4=[2*abs(fz4)/TD].^2;
    
    DIM=DIMs(dimv);
    dx4=[fx4(DIM,:);  geomean(fx4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz*0]],:))];
    dy4=[fy4(DIM,:);  geomean(fy4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz*0]],:))];
    dv4=[fv4(DIM,:);  geomean(fv4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz*0]],:))];
    dw4=[fw4(DIM,:);  geomean(fw4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz*0]],:))];
    du4=[fu4(DIM,:);  geomean(fu4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz*0]],:))];
    dz4=[fz4(DIM,:);  geomean(fz4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz*0]],:))];
    
    Bv=1e4;
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dx4'); px4=[sum(A<=0)+1]/[Bv+1];mx4=O;sx4=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dy4'); py4=[sum(A<=0)+1]/[Bv+1];my4=O;sy4=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dw4'); pw4=[sum(A<=0)+1]/[Bv+1];mw4=O;sw4=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dv4'); pv4=[sum(A<=0)+1]/[Bv+1];mv4=O;sv4=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,du4'); pu4=[sum(A<=0)+1]/[Bv+1];mu4=O;su4=std(A);
    [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dz4'); pz4=[sum(A<=0)+1]/[Bv+1];mz4=O;sz4=std(A);
    
    ps(:,cc,dimv)=[py4 px4 pw4 pv4 pu4 pz4];
    ms(:,cc,dimv)=[my4 mx4 mw4 mv4 mu4 mz4];
    ss(:,cc,dimv)=[sy4 sx4 sw4 sv4 su4 sz4];
  end
end
for cc=1:4
  [h, crit_p, ps(1,cc,:)]=fdr_bh(squeeze(ps(1,cc,:)));
  [h, crit_p, ps(2,cc,:)]=fdr_bh(squeeze(ps(2,cc,:)));
  [h, crit_p, ps(3,cc,:)]=fdr_bh(squeeze(ps(3,cc,:)));
  [h, crit_p, ps(4,cc,:)]=fdr_bh(squeeze(ps(4,cc,:)));
  [h, crit_p, ps(5,cc,:)]=fdr_bh(squeeze(ps(5,cc,:)));
  [h, crit_p, ps(6,cc,:)]=fdr_bh(squeeze(ps(6,cc,:)));
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
%Correlation across 6 neural/ocular measures (rows) in experiment 2  in 3 experimental conditions
figure
clrs={[0 0 0],[1 1 1],[0 0 0]};
msps={'o','^','+'};
for cc=1:3
  cdt=cdts{cc};
  clr=clrs{cc};
  msp=msps{cc};
  fy4=fft(mean(y4(wind,:,cdt),3));fy4=[2*abs(fy4)/TD].^2;
  fx4=fft(mean(x4(wind,:,cdt),3));fx4=[2*abs(fx4)/TD].^2;
  fw4=fft(mean(w4(wind,:,cdt),3));fw4=[2*abs(fw4)/TD].^2;
  fv4=fft(mean(v4(wind,:,cdt),3));fv4=[2*abs(fv4)/TD].^2;fv4=fv4+sqrt(mean(fv4(:).^2))/1000;
  
  DIM=11;
  dx4=fx4(DIM,:)./geomean(fx4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  dy4=fy4(DIM,:)./geomean(fy4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  dv4=fv4(DIM,:)./geomean(fv4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  dw4=fw4(DIM,:)./geomean(fw4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  dx4s(:,cc)=dx4;
  dy4s(:,cc)=dy4;
  dv4s(:,cc)=dv4;
  dw4s(:,cc)=dw4;
  
  subplot 241;hold on
  plot(dx4(1,:),dv4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('vEOG');ylabel('blink');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);

  subplot 242;hold on
  plot(dw4(1,:),dv4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('hEOG');ylabel('blink');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);

  subplot 243;hold on
  plot(dx4(1,:),dy4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('vEOG');ylabel('EEG');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);

  subplot 244;hold on
  plot(dw4(1,:),dy4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('hEOG');ylabel('EEG');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);
  
  [c(cc,1),pc(cc,1)]=corr(dx4(1,:)',dv4(1,:)');
  [c(cc,2),pc(cc,2)]=corr(dw4(1,:)',dv4(1,:)');
  [c(cc,3),pc(cc,3)]=corr(dx4(1,:)',dy4(1,:)');
  [c(cc,4),pc(cc,4)]=corr(dw4(1,:)',dy4(1,:)');
end

clrs={[0 0 0],[1 1 1],[0 0 0]};
msps={'o','^','+'};
pv1=finv(1-0.05,2,4*ffz);
pv2=finv(1-0.05,2,4*ffz);
for cc=1:3
  cdt=cdts{cc};
  clr=clrs{cc};
  msp=msps{cc};
  fw4=fft(mean(w4(wind,:,cdt),3));fw4=[2*abs(fw4)/TD].^2;
  fx4=fft(mean(x4(wind,:,cdt),3));fx4=[2*abs(fx4)/TD].^2;
  fu4=fft(mean(u4(wind,:,cdt),3));fu4=[2*abs(fu4)/TD].^2;fu4=fu4+sqrt(mean(fu4(:).^2))/1000;
  fz4=fft(mean(z4(wind,:,cdt),3));fz4=[2*abs(fz4)/TD].^2;
  
  DIM=11;
  dx4=fx4(DIM,:)./geomean(fx4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  dw4=fw4(DIM,:)./geomean(fw4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  dz4=fz4(DIM,:)./geomean(fz4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  du4=fu4(DIM,:)./geomean(fu4([[DIM-ffz]:[DIM-1],[DIM+1]:[DIM+ffz]],:));
  dx4s(:,cc)=dx4;
  dw4s(:,cc)=dw4;
  dz4s(:,cc)=dz4;
  du4s(:,cc)=du4;
  
  subplot 247;hold on
  plot(dx4(1,:),dz4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('vEOG');ylabel('pupil');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);

  subplot 248;hold on
  plot(dw4(1,:),dz4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('hEOG');ylabel('pupil');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);

  subplot 245;hold on
  plot(dx4(1,:),du4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('vEOG');ylabel('saccade');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);

  subplot 246;hold on
  plot(dw4(1,:),du4(1,:),msp,'color','k','markerfacecolor',clr)
  set(gca,'yscale','log');set(gca,'xscale','log');
  axis([0.01 3e2 0.01 3e2]);set(gca,'xtick',[0.1 1 10 1e2]);set(gca,'ytick',[0.1 1 1e2])
  xlabel('hEOG');ylabel('EEG');
  set(gca,'xtick',[.1 1 10 100]);set(gca,'ytick',[.1 1 10 100]);
  
  try
  [c(cc,1+4),pc(cc,1+4)]=corr(dx4(1,:)',dz4(1,:)');
  [c(cc,2+4),pc(cc,2+4)]=corr(dw4(1,:)',dz4(1,:)');
  [c(cc,3+4),pc(cc,3+4)]=corr(dx4(1,:)',du4(1,:)');
  [c(cc,4+4),pc(cc,4+4)]=corr(dw4(1,:)',du4(1,:)');
  end
end

[h, crit_p, adj_pc]=fdr_bh(pc);
adj_pc
c