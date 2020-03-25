function  PMOEA_MOEA_P(path,name,real_path,overlapping,isRealWorld,startIndex,endIndex,c)
global total_algo_time time_MOEA

global edgeMatrix 

global avgDegree

% Problem='社团检测';

%%%%%
name1=name;
tic;
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
             AdjMatrix = reIndex(AdjMatrix);        %%修改时后加的，编号重排YHP
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
         save([networkDataRoot, strcat(name1,'.mat')], 'AdjMatrix', '-v7.3');
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
t_0 = clock;





clear EDGE EDGE_Data edgeMatrix1 edgeMatrix2;
% Nodes = single(1:length(AdjMatrix));
% nodesIndex = true(1,length(AdjMatrix));

AdjMatrix = sparse(logical(AdjMatrix));
degree = single(full(sum(AdjMatrix)));
% degree = single(sum(AdjMatrix));

% if isempty(gcp('nocreate'))
%     parpool(4);            %% 开启并行池——YHP
% end

% sum_adj = sum(sum(AdjMatrix));
% avgDegree = floor(sum(degree)/length(AdjMatrix));

%% 稿件修改时加入：应对孤立点等异常情况的预处理——YHP

index=find(degree==0);
AdjMatrix(index,:)=[];
AdjMatrix(:,index)=[];
Nodes = single(1:length(AdjMatrix));
nodesIndex = true(1,length(AdjMatrix));
numVar = size(AdjMatrix, 1);

for i=1:numVar
    AdjMatrix(i, i) = 0;
end

degree = single(full(sum(AdjMatrix)));
sum_adj = full(sum(sum(AdjMatrix)));
avgDegree = floor(sum(degree)/length(AdjMatrix));

%%
 
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
    coreNodes = single(sort(findCoreNodesByDegree(AdjMatrix,degree,avgDegree))); 
    save([networkDataRoot, 'coreNodes'], 'coreNodes', '-v7.3');
    
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

% community = cell(1,length(currentCoreNodes));
%% 对每一个核心点进行进化扩充
networkDataRoot = sprintf('adjMatrix_coreNodes/%s/coreCommunity_%s/',name1,num2str(c));
if ~isdir(networkDataRoot) %判断路径是否存在
    mkdir(networkDataRoot);
end
cellCommunity =  cell(1,length(currentCoreNodes));



coreNode_number = length(currentCoreNodes);
save_root = ['result_statistics/',name,'/coreNode_number'];
  if ~isdir(save_root) %判断路径是否存在
             mkdir(save_root);    
    end
dlmwrite([save_root,'/coreNode_number_',num2str(c),'.txt'],[coreNode_number currentCoreNodes],' ');

% return;

t_1 = clock;
if isempty(gcp('nocreate'))
    parpool(50);            %% 开启并行池——YHP
end
t_p = clock;

parfor ci = 1:length(currentCoreNodes) %% 所有MOEA并行运行
%for ci = 1:length(currentCoreNodes) %% 所有MOEA并行运行
 [subPop_ci,idealp] = initialPopulation(currentCoreNodes(1,ci),AdjMatrix,sum_adj);
                        coreNode = currentCoreNodes(1,ci);
                        [weights,neighbors] = init_weight(popsize, niche);
                        chromosomes = subPop_ci;
                        for Gene = 1:max_gen
                            for i=1:popsize
                                i_neighbor_index = neighbors(i,:);
                                i_neighbor_chromosome = chromosomes(i_neighbor_index);
                                %% 交叉变异产生新解
                                child = crossover_mutation(coreNode,i_neighbor_chromosome,crossover_posibility,mutation_posibility,AdjMatrix,Nodes,nodesIndex);           
                                %% 评价新解 
                                child = evaluate(coreNode,child,sum_adj,AdjMatrix,coreNodes,Nodes,degree);
                                %% 更新参考点
                                 for h=1:2
                                     if child{1}(end-(2-h))<idealp(h) 
                                        idealp(h)=child{1}(end-(2-h));       
                                     end
                                 end
                                %% 更新邻居域
                               chromosomes=update_neighbour(idealp,chromosomes,child,i_neighbor_index,weights,niche); 
                            end
%                                 clc;
%                                 fprintf('%s第%2s轮,%5s问题,第%2s/%2s维,已完成%4s%%,耗时%5s秒\n',name,num2str(1),Problem,num2str(ci),num2str(length(coreNodes)),num2str(roundn(Gene/max_gen*100,-1)),num2str(roundn(toc,-2)));
                        end
%                          coreCommunityTime = [coreCommunityTime;roundn(toc,-2)];
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
%  community{ci} = chromosomes; 
 cellCommunity{ci} = chromosomes;

%  networkDataRoot3 = sprintf('results/%s/coreCommunityAll/coreCommunityAll_%s/',name1,num2str(c));
% if ~isdir(networkDataRoot3) %判断路径是否存在
%     mkdir(networkDataRoot3);
% end

% pathTime=sprintf('results/%s/coreCommunityTime/coreCommunityTime(%s)_%s.txt',name1,num2str(length(coreCommunityTime)),num2str(c));
% coreCommunityTime = [coreCommunityTime(1,1);coreCommunityTime(1,2:end)-coreCommunityTime(1,1:end-1)];
% savedata1(pathTime,[coreCommunityTime;0;mean(coreCommunityTime)]);

% coreCommunityName = sprintf('%s_%s_%s.mat',num2str(startIndex+ci-1),num2str(coreNode),name1);
% save([networkDataRoot2, coreCommunityName], 'chromosomes');

% coreCommunityName3 = sprintf('%s_%s_%s.mat',num2str(startIndex+ci-1),num2str(coreNode),name1);
% save([networkDataRoot3, coreCommunityName3], 'chromosomesAll');
end
% pathTime=sprintf('results/%s/time/coreCommunityTime/coreCommunityTime(coreNodeNum=%s)_%s.txt',name1,num2str(length(coreCommunityTime)),num2str(c));
% coreCommunityTime = [coreCommunityTime(1,1);coreCommunityTime(2:end,1)-coreCommunityTime(1:end-1,1)];
% savedata1(pathTime,[coreCommunityTime;0;mean(coreCommunityTime)]);








 networkDataRoot2 = sprintf('results/%s/coreCommunityFirstFront/coreCommunity_%s/',name1,num2str(c));

 time_MOEA = toc
 total_algo_time = total_algo_time + time_MOEA;
 
 if ~isdir(networkDataRoot2) %判断路径是否存在
    mkdir(networkDataRoot2);
end
coreCommunityName = sprintf('%s.mat',name1);
save([networkDataRoot2, coreCommunityName], 'cellCommunity', '-v7.3');



t_2 = clock;


first_used_time = etime(t_1,t_0);
second_used_time = etime(t_2,t_1);
partime = etime(t_p,t_1);

save_root = ['result_statistics/',name,'/first_used_time'];
  if ~isdir(save_root) %判断路径是否存在
             mkdir(save_root);    
  end
dlmwrite([save_root,'/first_used_time_',num2str(c),'.txt'],first_used_time,' ');


save_root = ['result_statistics/',name,'/second_used_time'];
  if ~isdir(save_root) %判断路径是否存在
             mkdir(save_root);    
  end
dlmwrite([save_root,'/second_used_time_',num2str(c),'.txt'],second_used_time,' ');



save_root = ['result_statistics/',name,'/paratera_time'];
  if ~isdir(save_root) %判断路径是否存在
             mkdir(save_root);    
  end
dlmwrite([save_root,'/paratera_time',num2str(c),'.txt'],partime,' ');

end

function [cells,idealp] = initialPopulation(coreNode,adjMatrix,sum_adj)
idealp = Inf*ones(1,2);
cells = cell(1,100);
for i  = 1:100
    adjIndex = single(find(adjMatrix(coreNode,:)==1));
    adjIndex =  single(adjIndex(randperm(length(adjIndex))));
    selectAdjNode =  single(adjIndex(1,1:randi(length(adjIndex))));
%     selectAdjNode = adjIndex;
    notAdjNode =  single(setdiff(1:length(adjMatrix),[selectAdjNode coreNode]));
    nodes = [coreNode selectAdjNode];
    d_in = full(sum(adjMatrix(nodes,nodes)));
    d_out = full(sum(adjMatrix(notAdjNode,nodes)));
    index = find(d_out ~= 0);
    p = rand(1,sum(d_out ~= 0));
    add_p2 = d_in(d_out ~= 0)./d_out(d_out ~= 0);
    nodes(index(p>add_p2)) = [];
    d_in = sum(d_in);%% 内部边数（/2）?团内节点的度(不/2)
    d_out = sum(d_out);%% 团间边数
    din_dout = d_in+d_out;
    ce = [unique([nodes coreNode]) 1 d_out/min(din_dout,sum_adj-din_dout)];
%   cell = [unique([nodes coreNode]) 1 -sum(d_in)/(sum(d_in)+sum(d_out)) ];
    for h=1:2
        idealp(h)=ce(end-(2-h));       %更新参考点---3.7
    end
    cells(1,i) = {ce};
end
end






