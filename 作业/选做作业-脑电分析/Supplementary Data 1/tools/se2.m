function s=se2(x,dim)
% function s=se2(x,dim)
% standard error
s=std(x,[],dim)/sqrt(size(x,1));