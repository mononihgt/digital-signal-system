function y=filterB(h,x)
y=filter(h,1,x);
y=y(floor(length(h)/2+1):end,:);
