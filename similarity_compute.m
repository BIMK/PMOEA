function [similarity] = similarity_compute(adjMatrix,edgeMatrix,edgeAdjMatrix)
%adjMatrix原始图的顶点邻接矩阵
%edgeMatrix边表
%edgeAdjMatrix图的边邻接矩阵
     similarity = zeros(size(edgeAdjMatrix));
    %[a,b] = find(edgeAdjMatrix==1);
    for it = 1:size(edgeAdjMatrix,1)
        for jt = 1:size(edgeAdjMatrix,2)
            if (edgeAdjMatrix(it,jt)==1)
                nodes1 = edgeMatrix(it,:);
                nodes2 = edgeMatrix(jt,:);
                nodes3 = setxor(nodes1,nodes2);
                %分别求nodes3集合中各个顶点的邻居（包括自身）
                adj_node1 = [find(adjMatrix(nodes3(1),:)),nodes3(1)];
                adj_node2 = [find(adjMatrix(nodes3(2),:)),nodes3(2)];
                intersect_adj_nodes = intersect(adj_node1,adj_node2);%邻居交集
                union_adj_nodes = union(adj_node1,adj_node2);        %邻居并集
                sim = size(intersect_adj_nodes,2)/size(union_adj_nodes,2);
                similarity(it,jt) = sim;

%                 if sim>=0.15                                    %%对相似度做阈值切割
%                     similarity(it,jt) = sim;
%                 else                                          %%同时对边的邻接矩阵也切割,相似度小，认为不邻接
%                     edgeAdjMatrix(it,jt) = 0;
%                 end
            end 
        end
    end
    
end