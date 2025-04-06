function md = meandf(X)
A1=angle(mean(X(:,1)));
A2=angle(mean(X(:,2)));
md=A1-A2;

end

