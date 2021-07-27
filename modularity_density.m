function [ y ] = modularity_density( clu_assignment,adj_mat,jj )
% ADJ_MAT: the adjacency matrix of the network.
% clu_assignment: the cluster label vector.
%%global ADJ_MAT
y = 0;
lamda = 0.5;
ADJ_MAT=adj_mat;% modularity density����Ĳ��������Li�����ġ�ȡ0.5��ʱ�����modularity density;
% ȡ1ʱ����ratio association��������С������; ȡ0ʱ����ratio cut�������ִ��������
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

