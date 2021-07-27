function  MOEAD_Net(path,name,real_path,overlapping,isRealWorld,startIndex,endIndex,c)
%%%%%
global edgeMatrix 
global AdjMatrix 
global coreNodes
global idealp 
global degree
global avgDegree
Problem='聚类问题';
%%%%%
name1=name;
tic
%% 自动根据读取的网络数据格式确定邻接矩阵和有无真是划分%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hasReal =  exist(real_path,'file');%% 有无真实划分 大于0，有
networkDataRoot = sprintf('adjMatrix_coreNodes/%s/',name1);
if ~isdir(networkDataRoot) %判断路径是否存在
    mkdir(networkDataRoot);
end
 networkData = sprintf('adjMatrix_coreNodes/%s/%s.mat',name1,name1);
  hasAdjMat = exist(networkData,'file');
 if hasAdjMat>0
     
    load(networkData);
 else
         AdjMatrix = single(load(path));
         AdjMatrix_size = size(AdjMatrix);
         if AdjMatrix_size(2)>2 %%邻接矩阵表示
        %     edgeNum = sum(sum(AdjMatrix));
            [edgeMatrix1,edgeMatrix2] = (find(AdjMatrix==1));
            edgeMatrix = [edgeMatrix2,edgeMatrix1];
         else                   %%边表表示
             edgeMatrix = AdjMatrix;
             needAddOne = 0;    %%是否需要加1
             numVar=(max(max(AdjMatrix(:,1)),max(AdjMatrix(:,2))));
             if find(AdjMatrix==0)>0 %% 从0开始编号
                  needAddOne = 1;
                  numVar=numVar+1;
             end
             edgeNum = AdjMatrix_size(1);
              AdjMatrix=Adjreverse(AdjMatrix,numVar,needAddOne);
         end
         save([networkDataRoot, strcat(name1,'.mat')], 'AdjMatrix');
 end

 numVar=single(size(AdjMatrix,1));
 if hasReal >0
     if overlapping == 0 %%非重叠
         Datalabel=(load(real_path));
         if size(Datalabel,2)==2  %%社团划分为“点--》社团”的2列形式
             Datalabel=(Datalabel(:,2)');
         end
     else                %%重叠
         if isRealWorld ==1 %%真实网络
              Datalabel=(load(real_path));
            realCommunity = label2community(load(real_path));
         else
            [realCommunity,~,~] = LFR_community2community(real_path);
            for k = 1:length(realCommunity)
                Datalabel(1,realCommunity{k}) = k;
            end
         end
     end
 else
     if overlapping == 0 %%非重叠
         Datalabel= false(1,numVar);
     else                %%重叠
         realCommunity = {};
     end
     
     
 end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear EDGE EDGE_Data edgeMatrix1 edgeMatrix2;

degree = single(sum(AdjMatrix));
sum_adj = sum(sum(AdjMatrix));
avgDegree = floor(sum(degree)/length(AdjMatrix));

 
strNetwork=name;
root = sprintf('results/%s/ParetoFront',name1,strNetwork);  %%建个文件保存实验数据
if ~isdir(root) %判断路径是否存在
    mkdir(root); 
end
root = sprintf('results/%s/metrics',name1,strNetwork);
if ~isdir(root) %判断路径是否存在
    mkdir(root);
end


M = 2;
popsize = 100;
niche =40;%邻居数量
max_gen=50;
crossover_posibility=0.9;%交叉概率
mutation_posibility=0.6;%变异概率


networkDataCoreNodes = sprintf('adjMatrix_coreNodes/%s/coreNodes.mat',name1);
hasCoreNodes =  exist(networkDataCoreNodes,'file');%% 有无保存的核心点数据
hasCoreNodes = 0;
if hasCoreNodes>0
        load(networkDataCoreNodes);
else
    coreNodes = single(sort(findCoreNodesByDegree(AdjMatrix))); 
    save([networkDataRoot, 'coreNodes'], 'coreNodes');
    
end
 
if endIndex<=0
    currentCoreNodes = coreNodes(1,startIndex:length(coreNodes));
else
    currentCoreNodes = coreNodes(1,startIndex:endIndex);   
end
  pathTime = sprintf('results/%s/time/coreCommunityTime/',name1,num2str(c));
    if ~isdir(pathTime) %判断路径是否存在
        mkdir(pathTime);
    end
    coreCommunityTime = [];


subPop = cell(length(currentCoreNodes),1);

community = {};
%% 对每一个核心点进行进化扩充
networkDataRoot = sprintf('adjMatrix_coreNodes/%s/coreCommunity_%s/',name1,num2str(c));
if ~isdir(networkDataRoot) %判断路径是否存在
    mkdir(networkDataRoot);
end
for ci = 1:length(currentCoreNodes)

 subPop{ci} = fff(currentCoreNodes(1,ci),AdjMatrix,sum_adj);
 idealp = Inf*ones(1,M);
                        coreNode = currentCoreNodes(1,ci);
                        [weights,neighbors] = init_weight(popsize, niche);
                        chromosomes = subPop{ci};
                        for Gene = 1:max_gen
                            for i=1:popsize
                                i_neighbor_index = neighbors(i,:);
                                i_neighbor_chromosome = chromosomes(i_neighbor_index);
                                %% 交叉变异产生新解
                                child = crossover_mutation(coreNode,i_neighbor_chromosome,crossover_posibility,mutation_posibility);           
%                                   child= crossover_mutation_TY(coreNode,i_neighbor_chromosome,crossover_posibility,mutation_posibility);
%                                   if(Gene<=30&&mod(Gene,10)==0)  
%                                      if rand<0.2
%                                          if ~isempty(child{1})
% %                                              random_len = randi(length(child{1}));
% %                                              for c_i = 1:randi(random_len)
%                                              for c_i = 1:length(child{1})
%                                                  subCommunity = find_k_complete(AdjMatrix,child{1}(1,c_i),3);
%                                                  child{1} = unique([child{1} subCommunity]);
%                                              end
%                                          end
%                                          
%                                      else
%                                             
%                                       integrateAdjNodes = setdiff(find(sum(AdjMatrix(child{1},:))>0),child{1});%%将上面所得的点看成一个整体，再以一定概率将这个整体的邻接点加入
%                                           integrateAdjNodes2 = integrateAdjNodes(randperm(length(integrateAdjNodes)));
%                                           if ~isempty(integrateAdjNodes2)
%                                                 child{1} = unique([child{1} integrateAdjNodes2(1,1:randi(length(integrateAdjNodes)))]);
%                                           end%3.4 %震荡pSim*avgSimilarity
%                                  
%                                      end
%                                  end
                               
                                   %% 评价新解 
                                child = evaluate(coreNode,child,sum_adj);
                                
                                %% 更新参考点
                                 for h=1:2
                                     if child{1}(end-(2-h))<idealp(h) 
                                        idealp(h)=child{1}(end-(2-h));       %更新参考点---3.7
                                     end
                                 end
                                %% 更新邻居域
                               chromosomes=update_neighbour(idealp,chromosomes,child,i_neighbor_index,weights,niche); %%更新种群%3.6
                            end
                               
                                clc;
                                fprintf('%s第%2s轮,%5s问题,第%2s/%2s维,已完成%4s%%,耗时%5s秒\n',name,num2str(1),Problem,num2str(ci),num2str(length(coreNodes)),num2str(roundn(Gene/max_gen*100,-1)),num2str(roundn(toc,-2)));

                        end
                      

                         coreCommunityTime = [coreCommunityTime;roundn(toc,-2)];
                        
                        %%种群个体去重
                        chromosomes = cellfun(@getArrayFromByteStream,cellfun(@uint8,containers.Map(cellfun(@char,cellfun(@getByteStreamFromArray,chromosomes,'un',0),'un',0),zeros(size(chromosomes))).keys,'un',0),'un',0);
                        %%非支配排序
                         objMat = zeros(length(chromosomes),2);
                         for i = 1:length(chromosomes)
                            objMat(i,:) = chromosomes{i}(1,end-1:end);
                         end
                         [FrontValue,~] = P_sort(objMat,'all');
                         chromosomesAll = chromosomes;
                        chromosomes = chromosomes(FrontValue==1);
 community{ci} = chromosomes; 
 networkDataRoot2 = sprintf('results/%s/coreCommunityFirstFront/coreCommunity_%s/',name1,num2str(c));
if ~isdir(networkDataRoot2) %判断路径是否存在
    mkdir(networkDataRoot2);
end
 networkDataRoot3 = sprintf('results/%s/coreCommunityAll/coreCommunityAll_%s/',name1,num2str(c));
if ~isdir(networkDataRoot3) %判断路径是否存在
    mkdir(networkDataRoot3);
end

% pathTime=sprintf('results/%s/coreCommunityTime/coreCommunityTime(%s)_%s.txt',name1,num2str(length(coreCommunityTime)),num2str(c));
% coreCommunityTime = [coreCommunityTime(1,1);coreCommunityTime(1,2:end)-coreCommunityTime(1,1:end-1)];
% savedata1(pathTime,[coreCommunityTime;0;mean(coreCommunityTime)]);

 coreCommunityName = sprintf('%s_%s.mat',name1,num2str(coreNode));
save([networkDataRoot2, coreCommunityName], 'chromosomes');

 coreCommunityName3 = sprintf('%s_%s.mat',name1,num2str(coreNode));
save([networkDataRoot3, coreCommunityName3], 'chromosomesAll');
end
pathTime=sprintf('results/%s/time/coreCommunityTime/coreCommunityTime(coreNodeNum=%s)_%s.txt',name1,num2str(length(coreCommunityTime)),num2str(c));
coreCommunityTime = [coreCommunityTime(1,1);coreCommunityTime(2:end,1)-coreCommunityTime(1:end-1,1)];
savedata1(pathTime,[coreCommunityTime;0;mean(coreCommunityTime)]);
end

function neighborSet = find_k_complete(adj,node,k)
% find one k-order complete subgraph in adj which contains node randomly

%     neighborSet = [];
    %---Modified by Tian, 7/27/2015---
    neighborSet = false(1,size(adj,1));
    %---------------------------------
    if length(find(adj(node,:))) > 1
        allSubMap = nchoosek(find(adj(node,:)),k-1);
        allSubMap = allSubMap(randperm(size(allSubMap,1)),:);
        for i = 1 : size(allSubMap,1)
            node = node(1,randperm(size(node,1)));
            nodes = [node,allSubMap(i,:)];
            if adj(nodes,nodes) + eye(k)
                 common1 = intersect(find_neighbors(adj,nodes(2)),find_neighbors(adj,nodes(3)));
                if length(common1)>1
%                 neighborSet = allSubMap(i,:);
%                 break;
                %---Modified by Tian, 7/27/2015---
                neighborSet(allSubMap(i,:)) = true;
                %---------------------------------
                end
            end
        end
    end
%     neighborSet = [neighborSet,node];
    %---Modified by Tian, 7/27/2015---
    neighborSet = [find(neighborSet),node];
    %---------------------------------
end

function cells = fff(coreNode,adjMatrix,sum_adj)
ce = {};
global idealp;
idealp = Inf*ones(1,2);
    adjIndex = single(find(adjMatrix(coreNode,:)==1));
    adjIndex =  single(adjIndex(randperm(length(adjIndex))));
%     selectAdjNode =  single(adjIndex(1,1:randi(length(adjIndex))));
    selectAdjNode = adjIndex;
    notAdjNode =  single(setdiff(1:length(adjMatrix),[selectAdjNode coreNode]));
    Nodes = [coreNode selectAdjNode];
    d_in = sum(adjMatrix(Nodes,Nodes));
    d_out = sum(adjMatrix(notAdjNode,Nodes));
    index = find(d_out ~= 0);
    p = rand(1,sum(d_out ~= 0));
    add_p2 = d_in(d_out ~= 0)./d_out(d_out ~= 0);
    Nodes(index(p>add_p2)) = [];
    d_in = sum(d_in);%% 内部边数（/2）?团内节点的度(不/2)
    d_out = sum(d_out);%% 团间边数
    din_dout = d_in+d_out;
    ce = [unique([Nodes coreNode]) 1 d_out/min(din_dout,sum_adj-din_dout)];
 
%   cell = [unique([Nodes coreNode]) 1 -sum(d_in)/(sum(d_in)+sum(d_out)) ];
    for h=1:2
        idealp(h)=ce(end-(2-h));       %更新参考点---3.7
    end
    cells = cell(1,100);
    cells(1,:) = {ce};
    
end


function neighbors = find_neighbors(adj,A)
% find all the neighbors of A in adj
    A = ismember(1:size(adj,1),A);
    neighbors = find(any(adj(A,:),1) & ~A);
end




