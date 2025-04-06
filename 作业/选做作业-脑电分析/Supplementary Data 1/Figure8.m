% **************************************************
% Figure 8
% **************************************************

clear;clc
fs=32;

figure;
Bv=1e4;
clr=[1 0 0;0 0 1;0 .7 0];
clf=[1 .5 .5;.5 .5 1;[.5 1 .5]];
for ind=1:3
  switch ind
    case 1
      load Data\TRF-yb;y=squeeze(y(:,5,:)); % time warpped responses, EEG
    case 2
      load Data\TRF-yeV; % time warpped responses, vEOG
    case 3
      load Data\TRF-yeH; % time warpped responses, hEOG
  end
  y=y/sqrt(mean(y(:).^2));
  TD=size(y,1);
  fy4=fft(y);fy4=[2*abs(fy4)/TD].^2;
  [O,Afy4]=bootstrap_for_vector(1e3,0.05,@geomean,fy4');
  
  f=1:size(fy4,1);f=f-1;f=f/size(fy4,1);f=f*fs;
  subplot(1,3,ind);hold on
  plot_shade(f(2:end),geomean(fy4(2:end,:)'),std(Afy4(:,2:end)),clf(ind,:))
  plot(f,geomean(fy4,2),'color',clr(ind,:),'linewidth',1.5)
  xlim([.7 4.25]);ylim([0 .35])
  set(gca,'ytick',0:.2:.8)
  
  DIM=10;
  dy4=[fy4(DIM,:); geomean(fy4([[DIM-2]:[DIM-1],[DIM+1]:[DIM+2]],:))];
  [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dy4'); p1=[sum(A<=0)+1]/[Bv+1];mx1=O;sx1=std(A);
  DIM=19;
  dy4=[fy4(DIM,:); geomean(fy4([[DIM-2]:[DIM-1],[DIM+1]:[DIM+2]],:))];
  [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dy4'); p2=[sum(A<=0)+1]/[Bv+1];mx2=O;sx2=std(A);
  DIM=37;
  dy4=[fy4(DIM,:); geomean(fy4([[DIM-2]:[DIM-1],[DIM+1]:[DIM+2]],:))];
  [O,A]=bootstrap_for_vector(Bv,0.05,@compare_geomean,dy4'); p4=[sum(A<=0)+1]/[Bv+1];mx4=O;sx4=std(A);
  
  ps(:,ind)=[p1 p2 p4];
  ms(:,ind)=[mx1 mx2 mx4];
  ss(:,ind)=[sx1 sx2 sx4];
end

%%
[h, crit_p, adj_p(1,:)]=fdr_bh(ps(1,:));
[h, crit_p, adj_p(2,:)]=fdr_bh(ps(2,:));
[h, crit_p, adj_p(3,:)]=fdr_bh(ps(3,:));
adj_p'

%%
for cc=1:3
  for dm=1:size(ms,1)
      txtv{cc,dm}=[num2str(adj_p(cc,dm)) ' (' num2str(ms(cc,dm)) '/' num2str(ss(cc,dm)) ')'];
    end
  end  
% end
txtv

%%
clear % ERP
Bv=1e3;
figure;
load Data\ERP-ye; % ERP, EOG
y=squeeze(ye(:,1,5:end,:));y=y(fs/4+[1:fs],:,:,:);y0=squeeze(ye(:,1,1,:));y0=y0(fs/4+[1:fs],:,:,:);
y=y-repmat(mean(mean(y(1:fs/4,:,:,:),2)),size(y,1),size(y,2),1,1);y=reshape(y,128,4,11,15);y=squeeze(mean(y,3));
t=1:size(y,1);t=t/fs;t=t-1/4;
subplot 132;hold on;
for sind=1:4
  [O,A]=bootstrap_for_vector(Bv,0.05,@mean,squeeze(y(:,sind,:))');
  plot_shade(t+0*(sind-1)/4,mean(y(:,sind,:),3)'-sind*2,std(A),min([.5 .5 1],1))
end
for sind=1:4
  plot(t+0*(sind-1)/4,mean(y(:,sind,:),3)-sind*2,'color',min([0 0 1],1),'linewidth',1.5)
end
for sind=1:4,plot(t+0*(sind-1)/4,mean(mean(y,2),3)-sind*2,':','color',[0 0 1],'linewidth',1.5);end
xlim([0 1]-1/4);set(gca,'xtick',0:1/4:1);ylim([-9.5 -.5]);grid
y2=y;

load Data\ERP-ye; % ERP, EOG
ye=ye/2;
y=squeeze(ye(:,2,5:end,:));y=y(fs/4+[1:fs],:,:,:);y0=squeeze(ye(:,2,1,:));y0=y0(fs/4+[1:fs],:,:,:);
y=y-repmat(mean(mean(y(1:fs/4,:,:,:),2)),size(y,1),size(y,2),1,1);y=reshape(y,128,4,11,15);y=squeeze(mean(y,3));
subplot 133;hold on
for sind=1:4
  [O,A]=bootstrap_for_vector(Bv,0.05,@mean,squeeze(y(:,sind,:))');
  plot_shade(t+0*(sind-1)/4,mean(y(:,sind,:),3)'-sind*2,std(A),[.5 1 .5])
end
for sind=1:4
  plot(t+0*(sind-1)/4,mean(y(:,sind,:),3)-sind*2,'color',[0 .7 0],'linewidth',1.5)
end
for sind=1:4,plot(t+0*(sind-1)/4,mean(mean(y,2),3)-sind*2,':','color',[0 .7 0],'linewidth',1.5);end
xlim([0 1]-1/4);set(gca,'xtick',0:1/4:1);ylim([-9.5 -.5]);grid
y3=y;

load Data\ERP-yb; % ERP, EEG
y=squeeze(yb(:,48,5:end,:));y=y(fs/4+[1:fs],:,:,:);y0=squeeze(yb(:,1,1,:));y0=y0(fs/4+[1:fs],:,:,:);
y=y-repmat(mean(mean(y(1:fs/4,:,:,:),2)),size(y,1),size(y,2),1,1);y=reshape(y,128,4,11,15);y=squeeze(mean(y,3));
subplot 131;hold on
for sind=1:4
  [O,A]=bootstrap_for_vector(Bv,0.05,@mean,squeeze(y(:,sind,:))');
  plot_shade(t+0*(sind-1)/4,mean(y(:,sind,:),3)'-sind*2,std(A),[1 .5 .5])
end
for sind=1:4
  plot(t+0*(sind-1)/4,mean(y(:,sind,:),3)-sind*2,'color',[1 0 0],'linewidth',1.5)
end
for sind=1:4,plot(t+0*(sind-1)/4,mean(mean(y,2),3)-sind*2,':','color',[1 0 0],'linewidth',1.5);end
xlim([0 1]-1/4);set(gca,'xtick',0:1/4:1);ylim([-9.5 -.5]);grid
y1=y;
