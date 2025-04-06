% **************************************************
% Figure 7
% **************************************************

clear;clc

addpath .\tools

fs=16;
load Data\yeE4
ye(abs(ye)>1e3)=0;
ye(:,1,:,:,:)=mean(ye(:,1:2,:,:,:),2);
ye(:,2,:,:,:)=mean(ye(:,3,:,:,:),2);
ye4=abs(ye);

load Data\yeE3
ye(abs(ye)>1e3)=0;
ye(:,1,:,:,:)=mean(ye(:,1:2,:,:,:),2);
ye(:,2,:,:,:)=mean(ye(:,3,:,:,:),2);
ye3=abs(ye);

load Data\yeE2
ye(abs(ye)>1e3)=0;
ye2=abs(ye);
ye2(:,:,:,:,1)=[ye2(:,:,:,:,1)+ye2(:,:,:,:,2)]/2;
ye2(:,:,:,:,2)=[];
ye2=ye2(:,:,:,:,[2 3 1 4]);

load Data\yeE1
ye(abs(ye)>1e3)=0;
ye1=abs(ye);

%%
figure
for eyech=1:2

    % Timecourse of the vertical EOG, and horizontal EOG power before, during, and after the stimulus. 
    t=0:1/fs:30;t=t-2;
    subplot(2,4,(eyech-1)*4+1);hold on
    
    plot_shade(t(1:size(ye1,1)),squeeze(mean(mean(abs(ye1(:,eyech,:,:,:).^2),3),4))',squeeze(se2(mean(abs(ye1(:,eyech,:,:,:).^2),3),4))',[.5 .5 .5]);
    plot(t(1:size(ye1,1)),squeeze(mean(mean(abs(ye1(:,eyech,:,:,:).^2),3),4)),'k');ylim([4e2 6e4]/(3*(eyech-1)+1))
    plot([0 0],[4e2 6e4]/(3*(eyech-1)+1),'y')
    plot([11 11],[4e2 6e4]/(3*(eyech-1)+1),'y')
    xlim([-1.5 t(size(ye1,1))-.5]);set(gca,'xtick',[0:2:10,11,12]);set(gca,'xgrid','on');set(gca,'yscale','log');
    if eyech==1
        ylabel('vEOG')
    elseif eyech==2
        ylabel('hEOG')
    end
    title('Ex1')
    Pre1(:,1)=mean(mean((ye1([t>-2 & t<0],eyech,:,:,:).^2),3),1);
    Mid1(:,1)=mean(mean((ye1([t>2 & t<11],eyech,:,:,:).^2),3),1);
    Pst1(:,1)=mean(mean((ye1([t>11 & t<14],eyech,:,:,:).^2),3),1);
    
    subplot(2,4,(eyech-1)*4+2);hold on
    plot_shade(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,1).^2),3),4))',squeeze(se2(mean(abs(ye2(:,eyech,:,:,1).^2),3),4))',[.5 .5 .5]);
    h21=plot(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,1).^2),3),4)),'k');ylim([4e2 6e4]/(3*(eyech-1)+1))
    plot([0 0],[4e2 6e4]/(3*(eyech-1)+1),'y')
    plot([12 12],[4e2 6e4]/(3*(eyech-1)+1),'y')
    xlim([-1.5 t(size(ye2,1))-.5]);set(gca,'xtick',[0:2:12]);set(gca,'xgrid','on');set(gca,'yscale','log')
    title('Ex2')
    Pre2(:,1:4)=mean(mean((ye2([t>-2 & t<0],eyech,:,:,:).^2),3),1);
    Mid2(:,1:4)=mean(mean((ye2([t>2 & t<12],eyech,:,:,:).^2),3),1);
    Pst2(:,1:4)=mean(mean((ye2([t>12 & t<14],eyech,:,:,:).^2),3),1);
    
    plot_shade(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,2).^2),3),4))',squeeze(se2(mean(abs(ye2(:,eyech,:,:,2).^2),3),4))',[.5 .5 1]);
    plot_shade(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,3).^2),3),4))',squeeze(se2(mean(abs(ye2(:,eyech,:,:,3).^2),3),4))',[.5 1 .5]);
    plot_shade(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,4).^2),3),4))',squeeze(se2(mean(abs(ye2(:,eyech,:,:,4).^2),3),4))',[1 .5 .5]);
    h22=plot(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,2).^2),3),4)),'--b');
    h23=plot(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,3).^2),3),4)),'--g');
    h24=plot(t(1:size(ye2,1)),squeeze(mean(mean(abs(ye2(:,eyech,:,:,4).^2),3),4)),'--r');ylim([4e2 6e4]/(3*(eyech-1)+1))
    xlim([-1.5 t(size(ye2,1))-.5]);set(gca,'xtick',[0:2:12]);set(gca,'xgrid','on');set(gca,'yscale','log')
    
    subplot(2,4,(eyech-1)*4+3);hold on
    plot_shade(t(1:size(ye3,1)),squeeze(mean(mean(abs(ye3(:,eyech,:,:,1).^2),3),4))',squeeze(se2(mean(abs(ye3(:,eyech,:,:,1).^2),3),4))',[.5 .5 .5]);
    plot_shade(t(1:size(ye3,1)),squeeze(mean(mean(abs(ye3(:,eyech,:,:,2).^2),3),4))',squeeze(se2(mean(abs(ye3(:,eyech,:,:,2).^2),3),4))',[1 .5 .5]);
    h31=plot(t(1:size(ye3,1)),squeeze(mean(mean(abs(ye3(:,eyech,:,:,1).^2),3),4)),'k');ylim([4e2 6e4]/(3*(eyech-1)+1))
    h32=plot(t(1:size(ye3,1)),squeeze(mean(mean(abs(ye3(:,eyech,:,:,2).^2),3),4)),'r');ylim([4e2 6e4]/(3*(eyech-1)+1))
    plot([0 0],[4e2 6e4]/(3*(eyech-1)+1),'y')
    plot([11 11],[4e2 6e4]/(3*(eyech-1)+1),'y')
    xlim([-1.5 t(size(ye3,1))-.5]);set(gca,'xtick',[0:2:10 11]);set(gca,'xgrid','on');set(gca,'yscale','log')
    title('Ex3')
    Pre3(:,1:2)=mean(mean((ye3([t>-2 & t<0],eyech,:,:,:).^2),3),1);
    Mid3(:,1:2)=mean(mean((ye3([t>2 & t<11],eyech,:,:,:).^2),3),1);
    Pst3(:,1:2)=mean(mean((ye3([t>11 & t<13],eyech,:,:,:).^2),3),1);
    
    subplot(2,4,(eyech-1)*4+4);hold on
    plot_shade(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,3).^2),3),4))',squeeze(se2(mean(abs(ye4(:,eyech,:,:,3).^2),3),4))',[.5 .5 .5]);ylim([4e2 2.5e4])
    plot_shade(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,4).^2),3),4))',squeeze(se2(mean(abs(ye4(:,eyech,:,:,4).^2),3),4))',[1 .5 .5]);ylim([4e2 2.5e4])
    h41=plot(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,3).^2),3),4)),'k');ylim([4e2 6e4]/(3*(eyech-1)+1))
    h42=plot(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,4).^2),3),4)),'r');ylim([4e2 6e4]/(3*(eyech-1)+1))
    xlim([-1.5 t(size(ye4,1))-.5]);set(gca,'xtick',[0:2:12]);set(gca,'xgrid','on');set(gca,'yscale','log')
    title('Ex4')
    
    plot_shade(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,1).^2),3),4))',squeeze(se2(mean(abs(ye4(:,eyech,:,:,1).^2),3),4))',[.5 .5 .5]);
    plot_shade(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,2).^2),3),4))',squeeze(se2(mean(abs(ye4(:,eyech,:,:,2).^2),3),4))',[1 .5 .5]);
    h43=plot(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,1).^2),3),4)),'--k');
    h44=plot(t(1:size(ye4,1)),squeeze(mean(mean(abs(ye4(:,eyech,:,:,2).^2),3),4)),'--r');ylim([4e2 6e4]/(3*(eyech-1)+1))
    plot([0 0],[4e2 6e4]/(3*(eyech-1)+1),'y')
    plot([12 12],[4e2 6e4]/(3*(eyech-1)+1),'y')
    xlim([-1.5 t(size(ye4,1))-.5]);set(gca,'xtick',[0:2:12]);set(gca,'xgrid','on');set(gca,'yscale','log')
    
    Pre4(:,1:4)=mean(mean((ye4([t>-2 & t<0],eyech,:,:,1:4).^2),3),1);
    Mid4(:,1:4)=mean(mean((ye4([t>2 & t<12],eyech,:,:,1:4).^2),3),1);
    Pst4(:,1:4)=mean(mean((ye4([t>12 & t<14],eyech,:,:,1:4).^2),3),1);
    
    % the significance of vertical EOG, and horizontal EOG power compared with the pre-stimulus and post-stimulus periods
    Bv=10000;
    [p_Pre1(:,:,eyech),p_Pst1(:,:,eyech)]=PreMidPst_bstrap(Pre1,Mid1,Pst1,Bv);
    [p_Pre2(:,:,eyech),p_Pst2(:,:,eyech)]=PreMidPst_bstrap(Pre2,Mid2,Pst2,Bv);  
    [p_Pre3(:,:,eyech),p_Pst3(:,:,eyech)]=PreMidPst_bstrap(Pre3,Mid3,Pst3,Bv);   
    [p_Pre4(:,:,eyech),p_Pst4(:,:,eyech)]=PreMidPst_bstrap(Pre4,Mid4,Pst4,Bv);  
    a1=mean([Pre1,Mid1,Pst1],1);
    a2=mean([Pre2,Mid2,Pst2],1);
    a3=mean([Pre3,Mid3,Pst3],1);
    a4=mean([Pre4,Mid4,Pst4],1);
    
end

%P-value, mean, and variance measures of vertical EOG, and horizontal EOG power
eyech=1;p_Pre11=p_Pre1(:,:,eyech),p_Pst11=p_Pst1(:,:,eyech),p_Pre21=p_Pre2(:,:,eyech),p_Pst21=p_Pst2(:,:,eyech)
p_Pre31=p_Pre3(:,:,eyech),p_Pst31=p_Pst3(:,:,eyech),p_Pre41=p_Pre4(:,:,eyech),p_Pst41=p_Pst4(:,:,eyech)

eyech=2;p_Pre12=p_Pre1(:,:,eyech),p_Pst12=p_Pst1(:,:,eyech),p_Pre22=p_Pre2(:,:,eyech),p_Pst22=p_Pst2(:,:,eyech)
p_Pre32=p_Pre3(:,:,eyech),p_Pst32=p_Pst3(:,:,eyech),p_Pre42=p_Pre4(:,:,eyech),p_Pst42=p_Pst4(:,:,eyech)

