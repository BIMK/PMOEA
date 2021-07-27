function child = evaluate(coreNode,child,sum_adj,AdjMatrix,coreNodes,Nodes,degree)
%     global AdjMatrix coreNodes Nodes degree;
%     notAdjNode = setdiff(1:length(AdjMatrix),child{1});
%     nodeIndex2 = false(1,length(Nodes));
%     nodeIndex2(1,child{1}) = 1;
%     notAdjNode = find((nodesIndex - nodeIndex2)>0);
    
    d_in = sum(single(full(sum(AdjMatrix(child{1},child{1})))));%% 内部边数（/2）?团内节点的度(不/2)
%     d_out = sum(single(sum(AdjMatrix(child{1},notAdjNode))));%% 团间边数
    d_out = sum(degree(1,child{1})) - d_in;
    
%     if ~isequal(d_out,d_out2)
%         stop = 1;
%     end
    
    
    sum_dindout = d_in+d_out;
%     containCurrentCoreNode = ismember(coreNode,child{1});
    currentNodesIndex = false(1,length(Nodes));
        currentNodesIndex(1,child{1}) = 1;
    containCurrentCoreNode =  currentNodesIndex(1,coreNode);
%       if containCurrentCoreNode22~=containCurrentCoreNode
%                 stop = 1;
%             end
    
    if containCurrentCoreNode == 1
        coreNodesIndex = false(1,length(Nodes));
        coreNodesIndex(1,coreNodes) = 1;
        
        
        coreNodeNum = sum(currentNodesIndex&coreNodesIndex);
        
%             if coreNodeNum~=num
%                 stop = 1;
%             end
        
        
        
    else
        coreNodeNum = inf;
    end
    child{1} = [child{1}  coreNodeNum  d_out/min(sum_dindout,sum_adj - sum_dindout)];
end
