% **************************************************
% Figure 9
% **************************************************

clear;clc
load Data\yeE1
ye(abs(ye)>1e3)=0;
b{1}=squeeze(mean(mean(ye(:,1,:,:,:),2),3));

load Data\yeE2
ye(abs(ye)>1e3)=0;
b{2}=squeeze(mean(mean(ye(:,1,:,:,:),2),3));

load Data\yeE3
ye(abs(ye)>1e3)=0;
b{3}=squeeze(mean(mean(ye(:,1:2,:,:,:),2),3));

load Data\yeE4
ye(abs(ye)>1e3)=0;
b{4}=squeeze(mean(mean(ye(:,1:2,:,:,:),2),3));

fs=128/8;
t=0:1/fs:20;t=t-2;
hB=fir1(fs,[.75 1.25]/(fs/2));

%%
% figure
for expno=1:4
  yy=b{expno};
%   subplot(1,4,expno);hold on
  for cds=1:size(yy,3)
    y=filterB(hB,yy(:,:,cds));
    if expno==1 || expno==2
      ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=1 & t<10],:),fs,9,size(yy,2)),2));
    elseif expno==3
      ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=1 & t<11],:),fs,10,size(yy,2)),2));
    elseif expno==4
      if cds<=2
        ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=4 & t<11],:),fs,7,size(yy,2)),2));
      else
        ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=8 & t<11],:),fs,3,size(yy,2)),2));end
    end
%     plot(t(1:size(y,1)),mean(y,2))
  end
%   axis([0 14 -45 45])
end
ys{2}(:,:,1)=[ys{2}(:,:,1)+ys{2}(:,:,2)]/2;
ys{2}(:,:,2)=[];
ys{2}=ys{2}(:,:,[2 3 1]);
ys{5}(:,:,1)=ys{4}(:,:,3);
ys{5}(:,:,2)=ys{4}(:,:,4);
ys{4}(:,:,3:4)=[];

%%
Bv=1000;
clr='krbm';
theta=0:pi/100:2*pi;
r=.7;
p=nan(5,4);
pv=nan(5,4,3);
figure
for expno=1:5
  z=ys{expno};
  for cds=1:size(z,3)
    p(expno,cds)=anova1(zscore(z(:,:,cds))',[],'off');
    [Od,Ad(:,expno,cds)]=bootstrap_for_vector(Bv,0.05,@maxmean,z(:,:,cds)'); 
    Ad(:,expno,cds)=2*pi*Ad(:,expno,cds)/fs;
    subplot(2,5,expno);hold on
    plot(mean(z(:,:,cds)'),clr(cds),'linewidth',2);axis([0 16 -9 9])
    set(gca,'xtick',[0:4:16]);set(gca,'xticklabel',[0:4:16]/fs)
    [dum,pkv]=max(mean(z(:,:,cds)'));
    errorbar(pkv,mean(z(pkv,:,cds)'),se(z(pkv,:,cds)'),clr(cds),'linewidth',2)

    subplot(2,5,expno+5);hold on
    pv(expno,cds,:)=prcRGplot(Ad(:,expno,cds),95,clr(cds),1.1-cds/10);cla
  end
end

%%
pv(3:5,1,:)=pv(3:5,1,:)+2*pi;
pv(3,2,:)=pv(3,2,:)+2*pi;
pv(4,2,:)=pv(4,2,:)+2*pi;
pv(5,2,:)=pv(5,2,:)+2*pi;
pv(2,:,:)=pv(2,:,:)+2*pi;
pv(1,:,:)=pv(1,:,:)+2*pi;
figure(100);hold on
for expno=1:5
  for cds=1:4
    pvnow=squeeze(pv(expno,cds,:))'/pi/2;
    if expno==5
      errorbar(4.4+cds/5-.1,pvnow(2),pvnow(2)-pvnow(1),pvnow(3)-pvnow(2),clr(cds),'linewidth',2)      
      plot(4.4+cds/5-.1,pvnow(2),['o' clr(cds)],'linewidth',2)
    else
      errorbar(expno+cds/5-.1,pvnow(2),pvnow(2)-pvnow(1),pvnow(3)-pvnow(2),clr(cds),'linewidth',2)
      plot(expno+cds/5-.1,pvnow(2),['o' clr(cds)],'linewidth',2)
    end
  end
end

%%
clear;clc
load Data\etkE4b
ys=squeeze(mean(mean(ys(:,1,:,:,:),2),3));
c{6}=ys;

load Data\etkE5
ys=squeeze(mean(mean(ys(:,1,:,:,:),2),3));
for cds=1:2,ys0(:,:,cds)=resample(ys(:,:,cds),5,4);end;fs0=20;
ys0=ys0([.2*4*fs0+1]:end,:,:,:,:,:,:);clear ys
for cds=1:2,ys(:,:,cds)=resample(ys0(:,:,cds),4,5);end;
c{5}=ys;

load Data\etkE4
c{4}=squeeze(mean(mean(ys(:,1,:,:,:),2),3));

load Data\etkE2
c{2}=squeeze(mean(mean(ys(:,1,:,:,:),2),3));
c{2}(:,:,1)=[c{2}(:,:,1)+c{2}(:,:,2)]/2;
c{2}(:,:,2)=[];
c{2}=c{2}(:,:,[2 1]);
clear ys

fs=128/8;
t=0:1/fs:30;t=t;
hB=fir1(fs,[.75 1.25]/(fs/2));

%%
% figure
for expno=[2 4 5 6]
  yy=c{expno};
%   subplot(1,5,expno);hold on
  for cds=1:size(yy,3)
    bsl=repmat(mean(yy(:,:,cds)),size(yy,1),1,1);
    y=filterBL(hB,yy(:,:,cds));
%     y=yy(:,:,cds);
    if expno==1 || expno==2
      ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=1 & t<10],:),fs,9,size(yy,2)),2));
    elseif expno==3
      ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=1 & t<11],:),fs,10,size(yy,2)),2));
    elseif expno==4
      if cds<=2
        ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=4 & t<11],:),fs,7,size(yy,2)),2));
      else
        ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=8 & t<11],:),fs,3,size(yy,2)),2));end
    elseif expno==5
      ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=2 & t<20],:),fs*2,9,size(yy,2)),2));
    elseif expno==6
      ys{expno}(:,:,cds)=squeeze(mean(reshape(y([t>=8 & t<22],:),fs*2,7,size(yy,2)),2));
    end
%     plot(t(1:size(y,1)),mean(y,2))
  end
%   axis([0 14 0 .4])
end

%%
Bv=1000;
clr='krbm';
theta=0:pi/100:2*pi;
r=.7;
p=nan(6,4);
pv=nan(6,4,3);
figure
for expno=[2 4 5 6]
  z=ys{expno};
  for cds=1:size(z,3)
    p(expno,cds)=anova1(zscore(z(:,:,cds))',[],'off');
    [Od,Ad(:,expno,cds)]=bootstrap_for_vector(Bv,0.05,@maxmean,z(:,:,cds)'); 
    if expno>=5
      Ad(:,expno,cds)=2*pi*Ad(:,expno,cds)/fs/2;
      Od=2*pi*Od/fs/2;
    else
      Ad(:,expno,cds)=2*pi*Ad(:,expno,cds)/fs;
      Od=2*pi*Od/fs;
    end
    
    if expno==5
      subplot(2,5,[expno-1]:5);hold on;
    elseif expno==6
      subplot(2,5,[[expno-2]:5]+2);hold on;
    else subplot(2,5,expno-1);hold on
    end
    plot(mean(z(:,:,cds)'),clr(cds),'linewidth',2);axis([0 16 0 .2])
    [dum,pkv]=max(mean(z(:,:,cds)'));
    errorbar(pkv,mean(z(pkv,:,cds)'),se(z(pkv,:,cds)'),clr(cds),'linewidth',2)
    if expno==5,xlim([0 32]);ylim([0 .4]);end
    set(gca,'xtick',[0:4:32]);set(gca,'xticklabel',[0:4:32]/fs)
    if expno==6,xlim([0 32]);ylim([0 .4]);end
    set(gca,'xtick',[0:4:32]);set(gca,'xticklabel',[0:4:32]/fs)
    
    if expno==6
      subplot(2,5,10);hold on
      pv(expno,cds,:)=prcRGplot(Ad(:,expno,cds),95,clr(cds),1.1-cds/10);
    elseif expno==5
      subplot(2,5,9);hold on
      pv(expno,cds,:)=prcRGplot(Ad(:,expno,cds),95,clr(cds),1.1-cds/10);
    else
      subplot(2,5,8);hold on
      pv(expno,cds,:)=prcRGplot(Ad(:,expno,cds),95,clr(cds),1.1-cds/10);
    end
    cla
    if Od-pv(expno,cds,2)>pi, Od=Od-2*pi;end
    if Od-pv(expno,cds,2)<-pi, Od=Od+2*pi;end
    pv(expno,cds,2)=Od;
  end
end
% subplot(255);ylim([0 .5])
% pv(5,2,:)=pv(5,2,:)-2*pi;

%%
pv(2,2,:)=pv(2,2,:)+2*pi;
pv(6,2,:)=pv(6,2,:)+2*pi;
pv(3:5,1,:)=pv(3:5,1,:)+2*pi;
pv(4,2,:)=pv(4,2,:)+2*pi;
pv(5,2,:)=pv(5,2,:)+2*pi;
pv(5,1,:)=pv(5,1,:)-2*pi;
% pv(5,:,:)=pv(5,:,:)*2;
pv(2,:,:)=pv(2,:,:)+2*pi;
figure(100);hold on
for expno=1:6
  for cds=1:4
    if expno==1 && cds==1,continue;end
    if expno>=5
      pvnow=squeeze(pv(expno,cds,:))'/pi;
    else
      pvnow=squeeze(pv(expno,cds,:))'/pi/2;
    end
    errorbar(expno+cds/5-.05,pvnow(2),pvnow(2)-pvnow(1),pvnow(3)-pvnow(2),clr(cds),'linewidth',2)
    plot(expno+cds/5-.05,pvnow(2),['x' clr(cds)],'linewidth',2)
  end
end
plot([0 6],[0 0],'b')
plot([0 4.9],[0 0]+1,'b')
plot([2.9 4.9],[0 0]+.25,'k:')
plot([4.9 5.9],[0 0]+.2,'k:')
plot([2.9 4.9],[0 0]+.5,'r')
plot([2.9 4.9],[0 0]+.75,'r:')
plot([4.9 5.9],[0 0]+1,'r')
plot([4.9 5.9],[0 0]+1.2,'r:')
for ix=.5:1:6
  plot([0 0]+ix+.4,[-6 6],'g:');end
xlim([.5 6.5]+.4)
set(gca,'xtick',1:5)
ylim([-.3 2])
