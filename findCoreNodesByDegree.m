function coreNodes = findCoreNodesByDegree(AdjMatrix,degree,avgDegree)
     degree2 = degree;
%      %%与核心点直接相连接的不考虑
%      adj_core_index = (sum(AdjMatrix(coreNodes,:)==1)>0);
%      %%核心点本身不再考虑
%      adj_core_index(:,coreNodes) = true;
%     degree2(1,adj_core_index) = -1;
    %%度小于平均度的不考虑
    degree2(degree2<=(avgDegree)) = -1;
    candidateCoreNodeNum = sum(degree2>0);
    [~,degree_index] = sort(degree2,'descend');%%去除核心点之后的度排序
    degree_index = degree_index(1,1:candidateCoreNodeNum);
%     coresNode = degree_index(1);
%     degree_index(1) = [];
    coreNodes = [];
    while ~isempty(degree_index)
        coreNodes = [coreNodes degree_index(1)];
        degree_index(1) = [];
        if isempty(degree_index)
            break;
        end
        adj = degree_index(AdjMatrix(coreNodes(end),degree_index)>0);
            for i = 1:length(adj)
                aaa = find(adj(i) == degree_index);
                if ~isempty(aaa)
                    degree_index(aaa) = [];
                end
            end
%         if find(find(AdjMatrix(coresNode(end),degree_index)>0) == degree_index)
%             
%         end
    end
    
    
    
    
    
%     
%     [~,index_AA] = sort(top_degree_nodes_coreNodes_distance(:,end),'descend');
%     if ceil(e)>length(index_AA)
%         len = length(index_AA);
%     else
%         len =ceil(e);
%     end
%     moreCoreNodes = degree_index(1,index_AA(1:len,:)');
%     coreNodes = [coreNodes moreCoreNodes];
end

