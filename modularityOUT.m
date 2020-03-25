function Q = modularity(solution,AdjacentMatrix,degree,edge_num )
%% ��ȡ������Ȩ�����ĳ�ֻ��ֵ�ģ���
%% ģ��ȵ�����ͼ��㷽������������:
%% �򻯰���Q
%%solution�����ǩ
%%degree �����ڵ�Ķ�degree=single(sum(AdjMatrix,2));
%%AdjacentMatrix�ڽӾ���
%%edge_num�������sum(degree)/2;

m=max(solution);
edge_num = edge_num*2;
Q=0;
for i=1:m
    index = single(find(solution==i));
    if ~isempty(index)
%         Community_AdjMatrix=AdjacentMatrix(index,index);
%         Degree_Matrix=degree(index)*degree(index)'/edge_num;
        d = degree(index);
        e = d';
        Q=Q+sum(sum(AdjacentMatrix(index,index)-d*e/edge_num,1));
    end
end
% clear Community_AdjMatrix Degree_Matrix index solution;
clear index solution d e;
Q=Q/edge_num;


% % % % % m=max(solution);
% % % % % Q=0;
% % % % % for i=1:m
% % % % %     index = solution==i;
% % % % %     if ~isempty(find(index, 1))
% % % % %     Community_AdjMatrix=AdjacentMatrix(index,index);
% % % % %     Degree_Matrix=degree(index)*degree(index)'/(2*edge_num);
% % % % %     Q=single(Q+sum(sum(Community_AdjMatrix-Degree_Matrix,1)));
% % % % %     end
% % % % % end
% % % % % clear Community_AdjMatrix Degree_Matrix;
% % % % % Q=Q/(2*edge_num);




end
