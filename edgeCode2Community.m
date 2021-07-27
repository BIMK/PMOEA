function community = edgeCode2Community(edgeCode,edgeAdjMatrix,edgeMatrix,AdjMatrix,similarity,M)
community = {};
nodeNum =size(AdjMatrix,1);
f_node = zeros(1,nodeNum);
edgeNum = size(edgeMatrix,1);

[ff,index] = sort(edgeCode);
    label_index = [ff;index];%第一行为边标签值，第二行为该标签所对应边在边矩阵里的索引
    labels = unique(edgeCode); 
    labelsNum = size(labels,2);
    labels_nodeSet = zeros(labelsNum,nodeNum);%存放边标签所对应的点集合
    labels_advSim = zeros(labelsNum,2);%存放边标签所对应的边的平均相似度和边数
    % % 解码
    if(labelsNum==edgeNum)%一个边一个划分，解码成一个点一个划分
        f_node(1,:) = 1:nodeNum;
    else
        for label_i = 1:labelsNum
           edgeIndex = label_index(2,find(label_index(1,:) == labels(label_i)));
           %计算该集合边之间的平均相似度
           sumSim = 0;
%           nchoosek(edgeIndex(1,:),2);
           for sim_i = 1:size(edgeIndex,2)
               for sim_j = (sim_i+1):size(edgeIndex,2)
                   sumSim = sumSim+similarity(sim_i,sim_j);
               end
           end
           labels_advSim(label_i,:) = [sumSim/size(edgeIndex,2),size(edgeIndex,2)];
            
            for label_j = 1:size(label_index,2)
                
                if(labels(label_i)==label_index(1,label_j))
                   edge = edgeMatrix(label_index(2,label_j),:);%由边索引找到边所包含的顶点
%                    f_node(edge(1)) = labels(label_i);%将边的标签值转化为对应点的标签值
%                    f_node(edge(2)) = labels(label_i);
                    %labels_nodeSet(label_i,1) = labels(label_i);
                    
                    labels_nodeSet(label_i,edge) = edge;%第i行表示第i个边标签包含的节点
                end
            end
        end
        
             
   
    %%化成community胞元数组
%     for ci = 1:size(labels_nodeSet,1)
%           community{ci} = labels_nodeSet(ci,find(labels_nodeSet(ci,:)~=0));
%     end
        for labels_nodeSet_i = 1:size(labels_nodeSet,1)
            for labels_nodeSet_j = (labels_nodeSet_i+1):size(labels_nodeSet,1)
               intersect_node = intersect(labels_nodeSet(labels_nodeSet_i,:),labels_nodeSet(labels_nodeSet_j,:));
               union_node = union(labels_nodeSet(labels_nodeSet_i,:),labels_nodeSet(labels_nodeSet_j,:));
               intersect_node(intersect_node==0) = [];
               union_node(union_node==0) = [];
               for k = 1:length(intersect_node)
                   set1 = labels_nodeSet(labels_nodeSet_i,:);
                   set2 = labels_nodeSet(labels_nodeSet_i,:);
                   f_in1 = sum(sum(AdjMatrix(find(set1(1,:)~=0),find(set1(1,:)~=0))))/2;
                   f_out1 = sum(sum(AdjMatrix(find(set1(1,:)~=0),find(set1(1,:)==0))));
                   set2(1,find(set1(1,:)==intersect_node(1,k)))=0;%社团A除去候选重叠点后的点集合
                   f_in2 = sum(sum(AdjMatrix(find(set2(1,:)~=0),find(set2(1,:)~=0))))/2;
                   f_out2 = sum(sum(AdjMatrix(find(set2(1,:)~=0),find(set2(1,:)==0))));
                    Q1 = f_in1/(f_in1+f_out1);
                    Q2 = f_in2/(f_in2+f_out2);
                   if Q1 <= Q2%模块度变小，删除该候选重叠点
                      labels_nodeSet(labels_nodeSet_i,find(labels_nodeSet(labels_nodeSet_i,:)==intersect_node(1,k)))=0;
                   end

                   set11 = labels_nodeSet(labels_nodeSet_j,:);
                   set22 = labels_nodeSet(labels_nodeSet_j,:);
                   f_in11 = sum(sum(AdjMatrix(find(set11(1,:)~=0),find(set11(1,:)~=0))))/2;
                   f_out11 = sum(sum(AdjMatrix(find(set11(1,:)~=0),find(set11(1,:)==0))));
                   set22(1,find(set11(1,:)==intersect_node(1,k)))=0;%社团A除去候选重叠点后的点集合
                   f_in22 = sum(sum(AdjMatrix(find(set22(1,:)~=0),find(set22(1,:)~=0))))/2;
                   f_out22 = sum(sum(AdjMatrix(find(set22(1,:)~=0),find(set22(1,:)==0))));
                    Q11 = f_in11/(f_in11+f_out11);
                    Q22 = f_in22/(f_in22+f_out22);
                   if Q11 <= Q22%模块度变小，删除该候选重叠点
                      labels_nodeSet(labels_nodeSet_j,find(labels_nodeSet(labels_nodeSet_j,:)==intersect_node(1,k)))=0;
                   end
                   
               end
               
               if abs(labels_advSim(labels_nodeSet_i,1)-labels_advSim(labels_nodeSet_j,1))>=0.15%比较两个标签的平均相似度
                            %将相似度较小的标签点集合中将公共部分剔除
                            if labels_advSim(labels_nodeSet_i,1)>labels_advSim(labels_nodeSet_j,1)
                                    setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_j,:),intersect_node);
                                    setdiffSet(setdiffSet==0)=[];
                                    labels_nodeSet(labels_nodeSet_j,:) = 0;
                                    labels_nodeSet(labels_nodeSet_j,setdiffSet) = setdiffSet;%只保留差集
                            elseif labels_advSim(labels_nodeSet_j,1)>labels_advSim(labels_nodeSet_i,1)
                                    setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_i,:),intersect_node);
                                    setdiffSet(setdiffSet==0)=[];
                                    labels_nodeSet(labels_nodeSet_i,:) = 0;
                                    labels_nodeSet(labels_nodeSet_i,setdiffSet) = setdiffSet;%只保留差集
                            end
%                    elseif labels_advSim(labels_nodeSet_i,2)>labels_advSim(labels_nodeSet_j,2)%比较两个标签的节点数目
%                             setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_j,:),intersect_node);
%                             setdiffSet(setdiffSet==0)=[];
%                             labels_nodeSet(labels_nodeSet_j,:) = 0;
%                             labels_nodeSet(labels_nodeSet_j,setdiffSet) = setdiffSet;%只保留差集
%                   elseif labels_advSim(labels_nodeSet_i,2)<=labels_advSim(labels_nodeSet_j,2)%比较两个标签的节点数目
%                             setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_i,:),intersect_node);
%                             setdiffSet(setdiffSet==0)=[];
%                             labels_nodeSet(labels_nodeSet_i,:) = 0;
%                             labels_nodeSet(labels_nodeSet_i,setdiffSet) = setdiffSet;%只保留差集
%                  else 
%                            if rand(1)>=.5 
%                                setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_j,:),intersect_node);
%                                setdiffSet(setdiffSet==0)=[];
%                                labels_nodeSet(labels_nodeSet_j,:) = 0;
%                                labels_nodeSet(labels_nodeSet_j,setdiffSet) = setdiffSet;%只保留差集
%                            else
%                                setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_i,:),intersect_node);
%                                setdiffSet(setdiffSet==0)=[];
%                                labels_nodeSet(labels_nodeSet_i,:) = 0;
%                                labels_nodeSet(labels_nodeSet_i,setdiffSet) = setdiffSet;%只保留差集
%                            end
               end
               
               
               
            end
        end
        
     end
    
     labels_nodeSet(find(sum(abs(labels_nodeSet),2)==0),:)=[];%除去全零的行
        for labels_nodeSet_j = 1:size(labels_nodeSet,1)
%             if labels_nodeSet(labels_nodeSet_j,nodeLabel_i)~=0
                community{labels_nodeSet_j} = labels_nodeSet(labels_nodeSet_j,find(labels_nodeSet(labels_nodeSet_j,:)~=0));
%             end
        end
%     for nodeLabel_i = 1:nodeNum
%         for labels_nodeSet_j = 1:size(labels_nodeSet,1)
%             if labels_nodeSet(labels_nodeSet_j,nodeLabel_i)~=0
%                 f_node(1,nodeLabel_i) = labels_nodeSet_j;
%             end
%         end
%     end
%     node_code = f_node;


end 