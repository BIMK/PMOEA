function [similarity] = similarity_compute(adjMatrix,edgeMatrix,edgeAdjMatrix)
%adjMatrixԭʼͼ�Ķ����ڽӾ���
%edgeMatrix�߱�
%edgeAdjMatrixͼ�ı��ڽӾ���
     similarity = zeros(size(edgeAdjMatrix));
    %[a,b] = find(edgeAdjMatrix==1);
    for it = 1:size(edgeAdjMatrix,1)
        for jt = 1:size(edgeAdjMatrix,2)
            if (edgeAdjMatrix(it,jt)==1)
                nodes1 = edgeMatrix(it,:);
                nodes2 = edgeMatrix(jt,:);
                nodes3 = setxor(nodes1,nodes2);
                %�ֱ���nodes3�����и���������ھӣ���������
                adj_node1 = [find(adjMatrix(nodes3(1),:)),nodes3(1)];
                adj_node2 = [find(adjMatrix(nodes3(2),:)),nodes3(2)];
                intersect_adj_nodes = intersect(adj_node1,adj_node2);%�ھӽ���
                union_adj_nodes = union(adj_node1,adj_node2);        %�ھӲ���
                sim = size(intersect_adj_nodes,2)/size(union_adj_nodes,2);
                similarity(it,jt) = sim;

%                 if sim>=0.15                                    %%�����ƶ�����ֵ�и�
%                     similarity(it,jt) = sim;
%                 else                                          %%ͬʱ�Աߵ��ڽӾ���Ҳ�и�,���ƶ�С����Ϊ���ڽ�
%                     edgeAdjMatrix(it,jt) = 0;
%                 end
            end 
        end
    end
    
end