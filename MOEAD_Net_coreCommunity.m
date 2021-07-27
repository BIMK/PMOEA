function  MOEAD_Net_coreCommunity(name,startIndex,endIndex)
%%%%%
global AdjMatrix 
global coreNodes
global idealp 

name1=name;
tic
% root = sprintf('results1/%s/%s/ParetoFront',name1,strNetwork);  %%建个文件保存实验数据
% if ~isdir(root) %判断路径是否存在
%     mkdir(root); 
% end
% root = sprintf('results1/%s/%s/metrics',name1,strNetwork);
% if ~isdir(root) %判断路径是否存在
%     mkdir(root);
% end
coreNodes = [];
networkData = sprintf('adjMatrix_coreNodes/%s/%s.mat',name1,name1);
load(networkData);
% AdjMatrix = AdjMatrix.AdjMatrix;
networkDataCoreNOdes = sprintf('adjMatrix_coreNodes/%s/coreNodes.mat',name1);
load(networkDataCoreNOdes);
% coreNodes = coreNodes.coreNodes;

coreNodes = coreNodes(1,startIndex:endIndex);
coreNodes
M = 2;
popsize = 100;
niche =40;%邻居数量
max_gen=50;
crossover_posibility=0.9;%交叉概率
mutation_posibility=0.6;%变异概率

% coreNodes = single(sort(findCoreNodes(AdjMatrix,2))); 
% coreNodes = single(sort(findCoreNodesByDegree(AdjMatrix))); 
% save([root, 'coreNodes'], 'coreNodes');

coreNodes

subPop = cell(length(coreNodes),1);

community = {};
%% 对每一个核心点进行进化扩充
for ci = 1:length(coreNodes)
    idealp = [];
    subPop{ci} = fff(coreNodes(1,ci),AdjMatrix);
idealp = Inf*ones(1,M);
                        coreNode = coreNodes(1,ci);
                        [weights,neighbors] = init_weight(popsize, niche);
                        chromosomes = subPop{ci};
                       for Gene = 1:max_gen
                            for i=1:popsize
                                i_neighbor_index = neighbors(i,:);
                                i_neighbor_chromosome = chromosomes(i_neighbor_index);
                                %% 交叉变异产生新解
                                [child,~] = crossover_mutation(coreNode,i_neighbor_chromosome,crossover_posibility,mutation_posibility);                                
                                  if(Gene<0.8*20)  
                                     if rand<0.2
                                         if ~isempty(child{1})
                                             for c_i = 1:length(child{1})
                                                 subCommunity = find_k_complete(AdjMatrix,child{1}(1,c_i),3);
                                                 child{1} = unique([child{1} subCommunity]);
                                             end
                                         end
                                         
                                         
                                      integrateAdjNodes = setdiff(find(sum(AdjMatrix(child{1},:))>0),child{1});%%将上面所得的点看成一个整体，再以一定概率将这个整体的邻接点加入
                                          integrateAdjNodes2 = integrateAdjNodes(randperm(length(integrateAdjNodes)));
                                          if ~isempty(integrateAdjNodes2)
                                                child{1} = unique([child{1} integrateAdjNodes2(1,1:randi(length(integrateAdjNodes)))]);
                                          end%3.4 %震荡pSim*avgSimilarity
                                    
                                     
                                 
                                     
                                     end
                                 end
                               
                                   %% 评价新解 
                                child = evaluate(coreNode,child);
                                
                                %% 更新参考点
                                 for h=1:2
                                     if child{1}(end-(2-h))<idealp(h) 
                                        idealp(h)=child{1}(end-(2-h));       %更新参考点---3.7
                                     else
                                     end
                                 end
                                %% 更新邻居域
                               chromosomes=update_neighbour(idealp,chromosomes,child,i_neighbor_index,weights); %%更新种群%3.6
                            end
                         
                                Problem='聚类问题';
                                M=2;
                                clc;
                                fprintf('%s第%2s轮,%5s问题,第%2s/%2s维,已完成%4s%%,耗时%5s秒\n',name,num2str(1),Problem,num2str(ci),num2str(length(coreNodes)),num2str(roundn(Gene/max_gen*100,-1)),num2str(roundn(toc,-2)));

                        end
                     
                        %%种群个体去重
                        chromosomes = cellfun(@getArrayFromByteStream,cellfun(@uint8,containers.Map(cellfun(@char,cellfun(@getByteStreamFromArray,chromosomes,'un',0),'un',0),zeros(size(chromosomes))).keys,'un',0),'un',0);
                        %%非支配排序
                         objMat = zeros(length(chromosomes),2);
                         for i = 1:length(chromosomes)
                            objMat(i,:) = chromosomes{i}(1,end-1:end);
                         end
                         [FrontValue,~] = P_sort(objMat,'all');
                        chromosomes = chromosomes(FrontValue==1);
 community{ci} = chromosomes;                       
end

networkDataRoot = sprintf('adjMatrix_coreNodes/%s/coreCommunity/',name1);
if ~isdir(networkDataRoot) %判断路径是否存在
    mkdir(networkDataRoot);
end
coreCommunityName = sprintf('%s_%s_%s',name1,num2str(startIndex),num2str(endIndex));
save([networkDataRoot, coreCommunityName], 'community');

 end 

% 
function cell = init_cell(cell,coreNode,adjMatrix)
global idealp;
idealp = Inf*ones(1,2);
    adjIndex = single(find(adjMatrix(coreNode,:)==1));
    adjIndex =  single(adjIndex(randperm(length(adjIndex))));
    selectAdjNode =  single(adjIndex(1,1:randi(length(adjIndex))));
%        selectAdjNode = adjIndex;
    notAdjNode =  single(setdiff(1:length(adjMatrix),[selectAdjNode coreNode]));
    d_in = sum(adjMatrix([coreNode selectAdjNode],[coreNode selectAdjNode]));
    d_out = sum(adjMatrix(notAdjNode,[coreNode selectAdjNode]));
    Nodes = [coreNode selectAdjNode];
    Nodes(d_out>d_in) = [];
    
    d_in = sum(d_in);%% 内部边数（/2）?团内节点的度(不/2)
    d_out = sum(d_out);%% 团间边数
    cell = [unique([Nodes coreNode]) 1 d_out/min((d_in+d_out),sum(sum(adjMatrix))-(d_in+d_out)) ];
 
%   cell = [unique([Nodes coreNode]) 1 -sum(d_in)/(sum(d_in)+sum(d_out)) ];
    for h=1:2
              idealp(h)=cell(end-(2-h));       %更新参考点---3.7
    end
end


function cells = fff(coreNode,adjMatrix)
ce = {};
global idealp;
idealp = Inf*ones(1,2);
    adjIndex = single(find(adjMatrix(coreNode,:)==1));
    adjIndex =  single(adjIndex(randperm(length(adjIndex))));
    selectAdjNode =  single(adjIndex(1,1:randi(length(adjIndex))));
%        selectAdjNode = adjIndex;
    notAdjNode =  single(setdiff(1:length(adjMatrix),[selectAdjNode coreNode]));
    d_in = sum(adjMatrix([coreNode selectAdjNode],[coreNode selectAdjNode]));
    d_out = sum(adjMatrix(notAdjNode,[coreNode selectAdjNode]));
    Nodes = [coreNode selectAdjNode];
    Nodes(d_out>d_in) = [];
    
    d_in = sum(d_in);%% 内部边数（/2）?团内节点的度(不/2)
    d_out = sum(d_out);%% 团间边数
    ce = [unique([Nodes coreNode]) 1 d_out/min((d_in+d_out),sum(sum(adjMatrix))-(d_in+d_out)) ];
 
%   cell = [unique([Nodes coreNode]) 1 -sum(d_in)/(sum(d_in)+sum(d_out)) ];
    for h=1:2
              idealp(h)=ce(end-(2-h));       %更新参考点---3.7
    end
    cells = cell(1,100);
    cells(1,:) = {ce};
    
end


function [ParetoFront1,remove]=find_error(ParetoFront1,AdjMatrix,CLique,Q)
        M=[];
        V=length(CLique);
        numVar=size(AdjMatrix,1);
        for in=1:size(ParetoFront1,1)
            change_node{in}=[];
            lable= ParetoFront1(in,1:numVar);
            %% %%后处理：2种方式-最大的Q lable= ParetoFront1(in,1:numVar);
            for i=1:V
                if(length(CLique{i})>=3)
                    for j=1:length(CLique{i})
                        m=CLique{i}(j);
                        k=lable(m);
                        index_last=find(lable==k);
                        A=setdiff(index_last,m);
                        neighbors=find(AdjMatrix(m,:));
                        com_max=multi_label(neighbors,lable);
                        ParetoFront1(in,m)=com_max;
                        if com_max~=k
                            change_node{in}=[change_node{in} m];
                        end
                    end
                end
            end
            QQ(in) = modularity(ParetoFront1(in,1:numVar),AdjMatrix);
            if QQ(in)>Q(in)                                                       %有疑问
                M=[M in];
            end
        end
        error_node=[];
        for i=1:length(M)
            error_node=[error_node change_node{M(i)}];
        end
        remove=unique(error_node);
        
    end


    function Clique=find_merge(Q,C_num,chromosomes,CLique)
        f(:,1)=-Q;
        f(:,2)=-C_num;
        V=length(CLique);
        FrontValue = P_sort(f,'all');   %front层
        B=find(FrontValue==1);
        K=chromosomes(B,1:V);
        [a,~]=size(K);
        MM=[];
        for i=1:a
            A=decode(K(i,:));
            MM=[MM;A];
        end
        visited=zeros(1,V);
        times=1;
        for i=1:V
            if visited(i)==0
                distance_i=zeros(1,V);
                for j=i+1:V
                    if visited(j)==0
                        A=MM(:,i)-MM(:,j);                           %求取差异性
                        index_0=find(A==0);
                        distance_i(j)=length(index_0)/a;     %%函数逻辑错误！
                    end
                end
                F=find(distance_i>0.9);
                visited(F)=1;
                
                Clique{times}=cell2mat(CLique([F i]));
                if length(Clique{times})>0
                    times=times+1;
                end
            end
        end
    end



    function Clique=find_clique(Clique,erase_node,numVar)
        t=length(Clique);
        for i=1:t
            Clique{i}=setdiff(Clique{i},erase_node);
        end
        R=[];
        for i=1:t
            R=[R Clique{i}];
        end
        R=setdiff(1:numVar,R);
        for i=1:length(R)
            t=t+1;
            Clique{t}=R(i);
        end
        Clique=Clique(find(cell2mat(cellfun(@(S)length(S),Clique,'UniformOutput',false))~=0));  %%剔除长度为0的子团！
        [~,rank] = sort(cell2mat(cellfun(@(s)s(1),Clique,'UniformOutput',false)));                %对元胞进行排序
        Clique= Clique(rank);
    end

    function Matrix=find_Matrix(Matrix,Clique)
        numVar=size(Matrix,1);
        degree=sum(Matrix,1);
        D=1:numVar;
        for i=1:length(Clique)
            C=Clique{i};
            vertex_min=find(D==C(1));
            while length(C)>1
                j=2;
                vertex_max=find(D==C(2));
                if length(vertex_max)==0
                    CD_node=find(Clique{i}==C(j));
                    Clique{i}(CD_node)=[];
                    C(j)=[];
                else
                    Matrix(:,vertex_min)=Matrix(:,vertex_min)+Matrix(:,vertex_max);
                    Matrix(vertex_min,:)=Matrix(vertex_min,:)+Matrix(vertex_max,:);
                    Matrix(vertex_min,vertex_min)=Matrix(vertex_min,vertex_min)-Matrix(vertex_min,vertex_max)-Matrix(vertex_max,vertex_max);
                    Matrix(vertex_max,:)=[];
                    Matrix(:,vertex_max)=[];
                    index=find(D==C(2));
                    D(index)=[];
                    C(j)=[];
                end
            end
        end
    end
    function f=inherit(ParetoFront1,Clique)
        RANK=cell2mat(cellfun(@(S)S(randi(length(S))),Clique,'UniformOutput',false));%%%RANK为1：V长的数组，里面存储哪些团在具体的个体里面的标签。
 
        RANK=ParetoFront1(RANK);    %为每个团的标签！充分利用进化得到的粒子信息！
    for j=1:max(RANK)
          O=find(RANK==j);
        if length(O)>0
         O=[O O(1)];
         O(1)=[];
        f(find(RANK==j))=O;
        end
    end
        
        
        
    end
    
    
    %%找k团
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

function neighbors = find_neighbors(adj,A)
% find all the neighbors of A in adj
    A = ismember(1:size(adj,1),A);
    neighbors = find(any(adj(A,:),1) & ~A);
end

