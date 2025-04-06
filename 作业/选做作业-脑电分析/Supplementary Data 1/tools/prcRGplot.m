function px=prcRGplot(x,pt,tagv,r)
if ~exist('r'),r=1;end
Nsamp=size(x,1);
Ncent=round(Nsamp*pt/100);
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

theta=px(1):0.01:px(2);
plot([0 cos(mx)]*r,[0 sin(mx)]*r,tagv,'linewidth',2)
plot([0 cos(theta) 0]*r,[0 sin(theta) 0]*r,tagv,'linewidth',2)

theta=0:0.01:2*pi;
plot([cos(theta)],[sin(theta)],'k')

if mx<px(1),mx=mx+2*pi;end
if mx>px(2),mx=mx-2*pi;end

px=[px(1) mx px(2)];
