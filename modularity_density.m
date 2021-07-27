function [ y ] = modularity_density( clu_assignment,adj_mat,jj )
% ADJ_MAT: the adjacency matrix of the network.
% clu_assignment: the cluster label vector.
%%global ADJ_MAT
y = 0;
lamda = 0.5;
ADJ_MAT=adj_mat;% modularity density里面的参数，详见Li的论文。取0.5的时候等于modularity density;
% 取1时等于ratio association，倾向发现小的社区; 取0时等于ratio cut，倾向发现大的社区。
clu_num = max(clu_assignment);
for i = 1:clu_num
    s_index = find(clu_assignment == i);
    s = ADJ_MAT(s_index,s_index);
    s_cardinality = length(s_index);
    
    ksum = 0;
    for j = 1:s_cardinality
        k = sum(ADJ_MAT(s_index(j),:));
        ksum = ksum + k;       
    end
    kin = sum(sum(s));
    kout = ksum - kin;
    
    ys = (2*lamda*kin - 2*(1-lamda)*kout)/s_cardinality;
    y = y + ys;
end 

end

