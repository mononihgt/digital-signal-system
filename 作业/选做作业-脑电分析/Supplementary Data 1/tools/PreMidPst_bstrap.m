function [p_Pre,p_Pst]=PreMidPst_bstrap(Pre,Mid,Pst,Bv)
% Bv=1000;
for block_ind=1:size(Pre,2)
    [O_PreM,A_PreM]=bootstrap_for_vector(Bv,0.05,@compare_geomean,[Pre(:,block_ind),Mid(:,block_ind)]);
    p_Pre(block_ind,1)=[sum(A_PreM<=0)+1]/[Bv+1];p_Pre(block_ind,2)=O_PreM;p_Pre(block_ind,3)=std(A_PreM);

    [O_PstM,A_PstM]=bootstrap_for_vector(Bv,0.05,@compare_geomean,[Pst(:,block_ind),Mid(:,block_ind)]);
    p_Pst(block_ind,1)=[sum(A_PstM<=0)+1]/[Bv+1]; p_Pst(block_ind,2)=O_PstM;p_Pst(block_ind,3)=std(A_PstM);
end

[~,~,p_Pre(:,1)]=fdr_bh(p_Pre(:,1));
[~,~,p_Pst(:,1)]=fdr_bh(p_Pst(:,1));


