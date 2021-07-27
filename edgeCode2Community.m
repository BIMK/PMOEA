function community = edgeCode2Community(edgeCode,edgeAdjMatrix,edgeMatrix,AdjMatrix,similarity,M)
community = {};
nodeNum =size(AdjMatrix,1);
f_node = zeros(1,nodeNum);
edgeNum = size(edgeMatrix,1);

[ff,index] = sort(edgeCode);
    label_index = [ff;index];%��һ��Ϊ�߱�ǩֵ���ڶ���Ϊ�ñ�ǩ����Ӧ���ڱ߾����������
    labels = unique(edgeCode); 
    labelsNum = size(labels,2);
    labels_nodeSet = zeros(labelsNum,nodeNum);%��ű߱�ǩ����Ӧ�ĵ㼯��
    labels_advSim = zeros(labelsNum,2);%��ű߱�ǩ����Ӧ�ıߵ�ƽ�����ƶȺͱ���
    % % ����
    if(labelsNum==edgeNum)%һ����һ�����֣������һ����һ������
        f_node(1,:) = 1:nodeNum;
    else
        for label_i = 1:labelsNum
           edgeIndex = label_index(2,find(label_index(1,:) == labels(label_i)));
           %����ü��ϱ�֮���ƽ�����ƶ�
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
                   edge = edgeMatrix(label_index(2,label_j),:);%�ɱ������ҵ����������Ķ���
%                    f_node(edge(1)) = labels(label_i);%���ߵı�ǩֵת��Ϊ��Ӧ��ı�ǩֵ
%                    f_node(edge(2)) = labels(label_i);
                    %labels_nodeSet(label_i,1) = labels(label_i);
                    
                    labels_nodeSet(label_i,edge) = edge;%��i�б�ʾ��i���߱�ǩ�����Ľڵ�
                end
            end
        end
        
             
   
    %%����community��Ԫ����
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
                   set2(1,find(set1(1,:)==intersect_node(1,k)))=0;%����A��ȥ��ѡ�ص����ĵ㼯��
                   f_in2 = sum(sum(AdjMatrix(find(set2(1,:)~=0),find(set2(1,:)~=0))))/2;
                   f_out2 = sum(sum(AdjMatrix(find(set2(1,:)~=0),find(set2(1,:)==0))));
                    Q1 = f_in1/(f_in1+f_out1);
                    Q2 = f_in2/(f_in2+f_out2);
                   if Q1 <= Q2%ģ��ȱ�С��ɾ���ú�ѡ�ص���
                      labels_nodeSet(labels_nodeSet_i,find(labels_nodeSet(labels_nodeSet_i,:)==intersect_node(1,k)))=0;
                   end

                   set11 = labels_nodeSet(labels_nodeSet_j,:);
                   set22 = labels_nodeSet(labels_nodeSet_j,:);
                   f_in11 = sum(sum(AdjMatrix(find(set11(1,:)~=0),find(set11(1,:)~=0))))/2;
                   f_out11 = sum(sum(AdjMatrix(find(set11(1,:)~=0),find(set11(1,:)==0))));
                   set22(1,find(set11(1,:)==intersect_node(1,k)))=0;%����A��ȥ��ѡ�ص����ĵ㼯��
                   f_in22 = sum(sum(AdjMatrix(find(set22(1,:)~=0),find(set22(1,:)~=0))))/2;
                   f_out22 = sum(sum(AdjMatrix(find(set22(1,:)~=0),find(set22(1,:)==0))));
                    Q11 = f_in11/(f_in11+f_out11);
                    Q22 = f_in22/(f_in22+f_out22);
                   if Q11 <= Q22%ģ��ȱ�С��ɾ���ú�ѡ�ص���
                      labels_nodeSet(labels_nodeSet_j,find(labels_nodeSet(labels_nodeSet_j,:)==intersect_node(1,k)))=0;
                   end
                   
               end
               
               if abs(labels_advSim(labels_nodeSet_i,1)-labels_advSim(labels_nodeSet_j,1))>=0.15%�Ƚ�������ǩ��ƽ�����ƶ�
                            %�����ƶȽ�С�ı�ǩ�㼯���н����������޳�
                            if labels_advSim(labels_nodeSet_i,1)>labels_advSim(labels_nodeSet_j,1)
                                    setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_j,:),intersect_node);
                                    setdiffSet(setdiffSet==0)=[];
                                    labels_nodeSet(labels_nodeSet_j,:) = 0;
                                    labels_nodeSet(labels_nodeSet_j,setdiffSet) = setdiffSet;%ֻ�����
                            elseif labels_advSim(labels_nodeSet_j,1)>labels_advSim(labels_nodeSet_i,1)
                                    setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_i,:),intersect_node);
                                    setdiffSet(setdiffSet==0)=[];
                                    labels_nodeSet(labels_nodeSet_i,:) = 0;
                                    labels_nodeSet(labels_nodeSet_i,setdiffSet) = setdiffSet;%ֻ�����
                            end
%                    elseif labels_advSim(labels_nodeSet_i,2)>labels_advSim(labels_nodeSet_j,2)%�Ƚ�������ǩ�Ľڵ���Ŀ
%                             setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_j,:),intersect_node);
%                             setdiffSet(setdiffSet==0)=[];
%                             labels_nodeSet(labels_nodeSet_j,:) = 0;
%                             labels_nodeSet(labels_nodeSet_j,setdiffSet) = setdiffSet;%ֻ�����
%                   elseif labels_advSim(labels_nodeSet_i,2)<=labels_advSim(labels_nodeSet_j,2)%�Ƚ�������ǩ�Ľڵ���Ŀ
%                             setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_i,:),intersect_node);
%                             setdiffSet(setdiffSet==0)=[];
%                             labels_nodeSet(labels_nodeSet_i,:) = 0;
%                             labels_nodeSet(labels_nodeSet_i,setdiffSet) = setdiffSet;%ֻ�����
%                  else 
%                            if rand(1)>=.5 
%                                setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_j,:),intersect_node);
%                                setdiffSet(setdiffSet==0)=[];
%                                labels_nodeSet(labels_nodeSet_j,:) = 0;
%                                labels_nodeSet(labels_nodeSet_j,setdiffSet) = setdiffSet;%ֻ�����
%                            else
%                                setdiffSet = setdiff(labels_nodeSet(labels_nodeSet_i,:),intersect_node);
%                                setdiffSet(setdiffSet==0)=[];
%                                labels_nodeSet(labels_nodeSet_i,:) = 0;
%                                labels_nodeSet(labels_nodeSet_i,setdiffSet) = setdiffSet;%ֻ�����
%                            end
               end
               
               
               
            end
        end
        
     end
    
     labels_nodeSet(find(sum(abs(labels_nodeSet),2)==0),:)=[];%��ȥȫ�����
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