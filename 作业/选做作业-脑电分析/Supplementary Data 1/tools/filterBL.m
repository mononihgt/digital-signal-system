function y=filterBL(h,x)
bsl=repmat(mean(x),size(x,1),1);
x=x-bsl;
y=filter(h,1,x);
y=y(floor(length(h)/2+1):end,:);
y=y+bsl(1:size(y,1),:);