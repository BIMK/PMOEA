function [f]= initialize_variables(subPop,adjMatrix,M)
%%  V�߽ڵ���Ŀ
global idealp;
global AdjMatrix;%ԭʼͼ�ڽӾ���
% adjMatrix = AdjMatrix;
nodes = size(AdjMatrix,1);
f_node = zeros(N,nodes);
f = zeros(N,V);
%% 
K = M + V;
for i = 1 : N
    f(i,V + 1: K) = evaluate_objective(f_node(i,1:nodes),AdjMatrix,ll);%�Ը���ĵı�������������
end
 idealp = min(f(:,V + 1: K));%�ο���
end

