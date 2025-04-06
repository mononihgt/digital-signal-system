function px=prcRGplot(x,pt,tagv,r)
if ~exist('r'),r=1;end
Nsamp=size(x,1);
Ncent=round(Nsamp*pt(1)/100);
x=mod(x,2*pi);
x=sort(x);
x=[x;x+2*pi];
for ind=1:Nsamp
  xr(ind,:)=x(Ncent+ind-1,:)-x(ind,:);
end
[minval,id]=min(xr);
for idx=1:size(x,2)
  px(idx,1)=x(id(idx),idx);
  px(idx,2)=x(Ncent+id(idx)-1,idx);
end
px=mod(px+pi,2*pi)-pi;
mx=mphase(x);
for idx=1:size(x,2)
  if px(idx,1)>px(idx,2)
    px(idx,2)=px(idx,2)+2*pi;
  end
  if px(idx,1)>mx
    px(idx,:)=px(idx,:)-2*pi;
  end
  if px(idx,2)<mx
    px(idx,:)=px(idx,:)-2*pi;
  end
end
px1=px;

% xx=sort(x(px1(1)<x & x<mx));
% px2(1)=prctile(xx,80);
% xx=x(sort(px1(2)>x & x>mx));
% px2(2)=prctile(xx,80);

theta=px2(1):0.01:px2(2);
rr=.8;
plot([0 cos(px1(1))*rr cos(theta) cos(px1(2))*rr 0]*r,[0 sin(px1(1))*rr sin(theta) sin(px1(2))*rr 0]*r,tagv,'linewidth',2)

plot([0 cos(mx)]*r,[0 sin(mx)]*r,tagv,'linewidth',2)

theta=0:0.01:2*pi;
plot([cos(theta)],[sin(theta)],'k')
