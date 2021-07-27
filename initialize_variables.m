function [f]= initialize_variables(subPop,adjMatrix,M)
%%  V边节点数目
global idealp;
global AdjMatrix;%原始图邻接矩阵
% adjMatrix = AdjMatrix;
nodes = size(AdjMatrix,1);
f_node = zeros(N,nodes);
f = zeros(N,V);
%% 
K = M + V;
for i = 1 : N
    f(i,V + 1: K) = evaluate_objective(f_node(i,1:nodes),AdjMatrix,ll);%对个体的的编码结果进行评价
end
 idealp = min(f(:,V + 1: K));%参考点
end

