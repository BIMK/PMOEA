function child = crossover_mutation(coreNode,i_neighbor_chromosome,crossover_posibility,mutation_posibility,AdjMatrix,Nodes,nodesIndex)
% global AdjMatrix;
% global Nodes;
% global nodesIndex;
child = {};
flag1to2 = true;
                %% ���Ϊ�գ�1�Ÿ��彻�浽2�Ÿ���
     if rand(1) <= crossover_posibility%%һ�����ʽ���
        if length(i_neighbor_chromosome)>=2%%�ھ��������������2
            index1 = randi(length(i_neighbor_chromosome));
            i_neighbor_1 = i_neighbor_chromosome(index1);%%����ھ������1
            i_neighbor_1 = {i_neighbor_1{1}(1:end-2)};
            i_neighbor_chromosome(index1) = [];%%ɾ���Ѿ�ѡ��ĸ���
            index2 = randi(length(i_neighbor_chromosome));
            i_neighbor_2 = i_neighbor_chromosome(index2);%%�����ʣ���ھ���ѡ�����2
            i_neighbor_2 = {i_neighbor_2{1}(1:end-2)};
            if length(i_neighbor_1{1})>length(i_neighbor_2{1})
                 nodes = setdiff(i_neighbor_1{1},i_neighbor_2{1});%%����1����2�
                 if ~isempty(nodes)
                     flag1to2 = true;  
                 end
            else
%                 nodes = setdiff(i_neighbor_2{1},i_neighbor_1{1});%%����1����2�
                nodesIndex1 = false(1,length(Nodes));
                nodesIndex2 = false(1,length(Nodes));
                nodesIndex1(1,i_neighbor_2{1}) = 1;
                nodesIndex2(1,i_neighbor_1{1}) = 1;
                 nodes = find((nodesIndex1 - nodesIndex2)>0);
%                       if ~isequal(nodes2222222222222,nodes)
%                           stop = 1;
%                       end
                
                
                if ~isempty(nodes)
                     flag1to2 = false;
                 end
            end

            if ~isempty(nodes)
                nodes2 = nodes(randperm(length(nodes)));
                newNode = nodes2(1:randi(length(nodes)));%%������ѡ��һ���ֽڵ�

                  if flag1to2
%                         notAdjNode = setdiff(1:length(AdjMatrix),[newNode i_neighbor_2{1}]);
                        nodesIndex2 = false(1,length(Nodes));
                        nodesIndex2(1,[newNode i_neighbor_2{1}]) = 1;
                        notAdjNode = find((nodesIndex - nodesIndex2)>0);
%                       if ~isequal(notAdjNode2222222222222,notAdjNode)
%                           stop = 1;
%                       end
                        
                        AAA = full(AdjMatrix(i_neighbor_2{1},newNode));
                        BBB = full(AdjMatrix(notAdjNode,newNode));
                        if size(AAA,1)>1
%                             d_in = sum(AdjMatrix(i_neighbor_2{1},newNode));
                            d_in = sum(AAA);
                        else
%                             d_in = AdjMatrix(i_neighbor_2{1},newNode);
                            d_in = AAA;
                        end
                        
                         if size(BBB,1)>1
%                              d_out = sum(AdjMatrix(notAdjNode,newNode));
                             d_out = sum(BBB);
                         else
%                              d_out = AdjMatrix(notAdjNode,newNode);
                             d_out = BBB;
                         end
                        
                          index = find(d_out ~= 0);
    p = rand(1,sum(d_out ~= 0));
     add_p2 = d_in(d_out ~= 0)./d_out(d_out ~= 0);
%      add_p2 = mapminmax(add_p,0,1);
    if isempty(p)
        p = [];         %%Ӧ��ά����һ�¡���YHP
    end
    if isempty(add_p2)
        add_p2 = [];    %%Ӧ��ά����һ�¡���YHP
    end
    newNode(index(p>add_p2)) = [];
                    
%                          newNode(d_out>d_in) = [];
                        child{1} = [i_neighbor_2{1} newNode];
                  else
%                       notAdjNode = setdiff(Nodes,[newNode i_neighbor_1{1}]);
                      
                      nodesIndex2 = false(1,length(Nodes));
                      nodesIndex2(1,[newNode i_neighbor_1{1}]) = 1;
                      notAdjNode = find((nodesIndex - nodesIndex2)>0);
                      
%                       if ~isequal(notAdjNode22222,notAdjNode)
%                           stop = 1;
%                       end
                      
                      
                      
                      
                      
                      
                 
                      AAA = full(AdjMatrix(i_neighbor_1{1},newNode));
                        BBB = full(AdjMatrix(notAdjNode,newNode));
                        if size(AAA,1)>1
%                             d_in = sum(AdjMatrix(i_neighbor_1{1},newNode));
                            d_in = sum(AAA);
                        else
%                             d_in = AdjMatrix(i_neighbor_1{1},newNode);
                            d_in = AAA;
                        end
                        
                         if size(BBB,1)>1
                            d_out = sum(BBB);
                         else
                             d_out = BBB;
                         end           
%                       
%                       d_in = sum(AdjMatrix(i_neighbor_1{1},newNode));
%                       d_out = sum(AdjMatrix(notAdjNode,newNode));
                      
                     index = find(d_out ~= 0);
    p = rand(1,sum(d_out ~= 0));
     add_p2 = d_in(d_out ~= 0)./d_out(d_out ~= 0);
%      add_p2 = mapminmax(add_p,0,1);
    if isempty(p)
        p = [];         %%Ӧ��ά����һ�¡���YHP
    end
    if isempty(add_p2)
        add_p2 = [];    %%Ӧ��ά����һ�¡���YHP
    end
    newNode(index(p>add_p2)) = [];
                      
%                        newNode(d_out>d_in) = [];

                        child{1} = [i_neighbor_1{1} newNode];
                  end
            else                                          %%�Ϊ��(i_neighbor_1��i_neighbor_2����ȫ�������߱�������ϵ��)
                if length(i_neighbor_1) > length(i_neighbor_2)
                    moreNodes = setdiff(i_neighbor_1{1},intersect(i_neighbor_1{1},i_neighbor_2{1}));
                    moreNodes = moreNodes(randperm(length(moreNodes)));
                    moreNodes = moreNodes(randi(length(moreNodes)));
                    
                    %% ��֤����Ľڵ� k_in>k_out
                    
                    %%
                    child{1} = [i_neighbor_2{1} moreNodes];
                elseif length(i_neighbor_1) < length(i_neighbor_2)
                    moreNodes = setdiff(i_neighbor_2{1},intersect(i_neighbor_1{1},i_neighbor_2{1}));
                    moreNodes = moreNodes(randperm(length(moreNodes)));
                    moreNodes = moreNodes(randi(length(moreNodes)));
                     %% ��֤����Ľڵ� k_in>k_out
                    
                    %%
                    child{1} = [i_neighbor_1{1} moreNodes];
                else%%���ɾ��һЩ�ڵ�
                    length1 = length(i_neighbor_1{1});
                    if length1 > 2
                     child{1} = [i_neighbor_1{1}(1:randi(length1)) coreNode];%%������ĵ� ��ֹ��ɾ��
                    else 
                     child{1} = i_neighbor_1{1};
                    end
               end
            end
        else%%�ھ��������� ���䶯 
            child{1} = i_neighbor_chromosome{1}(1:end-2);
        end
     end

     
    if ~isempty(child)&&length(child{1})>1
        if  rand(1) <= mutation_posibility%%һ�����ʱ���
%               integrateAdjNodes = setdiff(find(sum(AdjMatrix(child{1},:))>0),child{1});%%���������õĵ㿴��һ�����壬����һ�����ʽ����������ڽӵ����
              AAAAA = full(any(AdjMatrix(child{1},:)));
%               AAAAA_index = find(AAAAA>0);
%               integrateAdjNodes = setdiff(find(AAAAA>0),child{1});%%���������õĵ㿴��һ�����壬����һ�����ʽ����������ڽӵ���� 
              
              %% ���߼�������� ���Ƿ����������Ч��
              deleteNodes = zeros(1,length(Nodes),'logical');
              deleteNodes(1,child{1}) = 1;
              integrateAdjNodes = find((AAAAA-deleteNodes)>0);
              
%              if  ~isequal(integrateAdjNodes2,integrateAdjNodes)
%                  stop = 1;
%              end
              
              integrateAdjNodes2 = integrateAdjNodes(randperm(length(integrateAdjNodes)));
              if isempty(integrateAdjNodes2 )
                  selectedNodes = [];
              else
                   selectedNodes = integrateAdjNodes2(1:randi(length(integrateAdjNodes2)));%%���ѡ��һ���ֽڵ�
%                    notAdjNode = setdiff(Nodes,[selectedNodes child{1}]);
                   
                   
                    allNodesIndex = ones(1,length(Nodes),'logical');
                    setdiffNodes = zeros(1,length(Nodes),'logical');
                    setdiffNodes(1,[selectedNodes child{1}]) = 1;
                    notAdjNode = find((allNodesIndex - setdiffNodes)>0);
                    
%                     isequal(notAdjNode222222,notAdjNode)
%                     if  ~isequal(notAdjNode222222,notAdjNode)
%                         stop = 1;
%                     end
                    
                    AAA = full(AdjMatrix(child{1},selectedNodes));
                    BBB = full(AdjMatrix([notAdjNode selectedNodes],selectedNodes));
                    if size(AAA,1) >1
                         d_in = sum(AAA);
                    else
                         d_in = AAA;
                    end
                    
                    if size(BBB,1) >1
                         d_out = sum(BBB);
                    else
                         d_out = BBB;
                    end
%                     d_in = sum(AdjMatrix(child{1},selectedNodes));
%                     d_out = sum(AdjMatrix([notAdjNode selectedNodes],selectedNodes));
                    
                   index = find(d_out ~= 0);
                    p = rand(1,sum(d_out ~= 0));
                    add_p2 = d_in(d_out ~= 0)./d_out(d_out ~= 0);
%      add_p2 = mapminmax(add_p,0,1);
                    if isempty(p)
                        p = [];         %%Ӧ��ά����һ�¡���YHP
                    end
                    if isempty(add_p2)
                        add_p2 = [];    %%Ӧ��ά����һ�¡���YHP
                    end
                    selectedNodes(index(p>add_p2)) = [];
%                     selectedNodes(d_out>=4*d_in) = [];
              end
              %%�����ɾ��
            if rand(1) <= 0.3%%0.7 %%����
%                 if ~isempty(child{1})
%                                              for c_i = 1:length(child{1})
%                                                  subCommunity = find_k_complete(AdjMatrix,child{1}(1,c_i),3);
%                                                  child{1} = unique([child{1} subCommunity]);
%                                              end
%                 end
                community = [child{1} selectedNodes];
                 child{1} = community;
            else
                child{1} = child{1}(randperm(length(child{1})));
                deletedNode = child{1}(1:randi(length(child{1})));%%���ѡ��һ���ֽڵ�
                 child{1} = setdiff(child{1},deletedNode);
                child{1} = [child{1} coreNode];
            end
        end
    end
    if isempty(child)
        child{1} = i_neighbor_chromosome{1}(1:end-2);
    end
    child{1} = unique(child{1});
end