function  MOEAD_Net_getFinalCommunity2(path,name,real_path,overlapping,isRealWorld,c)
%%%%%
global edgeMatrix 
global AdjMatrix 
% global coreNodes
% global degree

%%%%%
name1=name;
AdjMatrix = [];
tic
%% 导入邻接矩阵
%% 自动根据读取的网络数据格式确定邻接矩阵和有无真是划分%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hasReal =  exist(real_path,'file');%% 有无真实划分 大于0，有
 networkData = sprintf('adjMatrix_coreNodes/%s/%s.mat',name1,name1);
 hasNetWorkAdj =  exist(networkData,'file');
if hasNetWorkAdj
    load(networkData);
else
    AdjMatrix = single(load(path));
end
networkData2 = sprintf('adjMatrix_coreNodes/%s/coreNodes.mat',name1);
load(networkData2);
if isempty(AdjMatrix)
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
              AdjMatrix = Adjreverse(AdjMatrix,numVar,needAddOne);
         end
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
 
%  degree=sum(AdjMatrix,1);
%  index=find(degree==0);
%  AdjMatrix(index,:)=[];
%  AdjMatrix(:,index)=[];
 
% degree = single(sum(AdjMatrix));
% avgDegree = floor(sum(degree)/length(AdjMatrix));


community = {};
networkDataRoot = sprintf('results/%s/coreCommunityFirstFront/coreCommunity_%s/',name1,num2str(c));%%获取第c次扩充的核心点
    path = networkDataRoot;  
    fileExt = '*.mat';  
    files = dir(fullfile(path,fileExt));  
    for i=1:size(files,1)  
        fileName = strcat(path,files(i,1).name);  
        coreCommunity = load(fileName);
        community = [community {coreCommunity.chromosomes}];
    end; 
    
     pathTime = sprintf('results/%s/time/postProcessTime/',name1);
    if ~isdir(pathTime) %判断路径是否存在
        mkdir(pathTime);
    end
    pathResult = sprintf('results/%s/Q_NMI/',name1);
    if ~isdir(pathResult) %判断路径是否存在
        mkdir(pathResult);
    end
    
    clock1 = clock;
    
    
    
    
    
    [communityIndex,bestQ] = run_ga(name,community);

    findNodes = [];
    for ii = 1:length(communityIndex)
        myCommunity{ii} = community{ii}{communityIndex(ii)}(1,1:end-2);
        findNodes = unique([findNodes myCommunity{ii}]);
    end
    
    myLabel = zeros(1,length(Datalabel));
    for ii = 1:length(myCommunity)
        myLabel(myCommunity{ii}) = ii;
    end
    modnmi0 = [];
    modnmi0(1,1) = bestQ ;%%不考虑离散点直接计算
    if sum(Datalabel) ~= 0
      modnmi0(1,2) = nmi(Datalabel,myLabel(1,:));
    else
      modnmi0(1,2) = 0;
    end



remainNodes = setdiff(1:length(AdjMatrix),findNodes);


%% 对剩余的节点进行分派，选取Q最大的
% for i = 1:length(remainNodes)
%     i_Q = [];
%     for j = 1:length(myCommunity)
%         myLabel(1,remainNodes(1,i)) = j;
%         i_Q = [i_Q modularity(myLabel,AdjMatrix)];
%     end
%     [maxQ,index] = max(i_Q);
%      myLabel(1,remainNodes(1,i)) = index;
%      myCommunity{index} = unique([ myCommunity{index} remainNodes(1,i)]);
% end
%% 查找邻居所在团
% for i = 1:length(remainNodes)
%     i_Q = [];
%     for j = 1:length(myCommunity)
%         if ~isempty(find(AdjMatrix(remainNodes(1,i),myCommunity{j})>0)) %%有邻居所在团
%             myLabel(1,remainNodes(1,i)) = j;%%加入邻居所在团
%             i_Q = [i_Q;j,modularity(myLabel,AdjMatrix)];
%         end
%     end
%     if ~isempty(i_Q)
%          [maxQ,index] = max(i_Q(:,2));
%         myLabel(1,remainNodes(1,i)) = i_Q(index,1);
%         myCommunity{i_Q(index,1)} = unique([ myCommunity{i_Q(index,1)} remainNodes(1,i)]);
%     end
%    
% end
%%%%%%%%%加入最多邻居所在社团
% % for i = 1:length(remainNodes)
% %     for j = 1:length(myCommunity)
% %         n = find(AdjMatrix(remainNodes(i),:)>0);
% %         n_c_num = zeros(1,length(myCommunity));
% %         if ~isempty(n)
% %             for k = 1:length(n)
% %                 if ismember(n(k),myCommunity{j})
% %                     n_c_num(j) =  n_c_num(j)+1;
% %                 end
% %             end
% %         end
% %     end
% %     [~,remain_i_label] = max(n_c_num);
% %     myLabel(1,remainNodes(1,i)) = remain_i_label;%%加入邻居所在团
% % end


% tic
%     for i = 1:length(remainNodes)
%         i_nebor = find(AdjMatrix(remainNodes(i),:)>0);
%         i_nebor_label = myLabel(i_nebor);
%         i_nebor_label(i_nebor_label==0) =[];
%         if ~isempty(i_nebor_label)
%             X = unique(i_nebor_label);
%             [M,N]=hist(i_nebor_label,X);
%             [~,index] = max(M);
%             myLabel(remainNodes(i)) = N(index);
%         else
%            myLabel(remainNodes(i)) = randi(length(myCommunity)); 
%         end
%     end
% toc











% %% 对剩余的节点进行分派，每个点一个团  然后用局部调整策略
% findCommunityLen = length(myCommunity);
% for i = 1:length(remainNodes)
%      myCommunity{findCommunityLen+i} = remainNodes(1,i);
% end

topTree = {};
 finalCommunity = {};
     treeIndex = 1;
     topTree{treeIndex} = myCommunity;

         COR = zeros(length(myCommunity),length(myCommunity));
         treeIndex = 1;
         for top_j = 1:length(myCommunity)-1
             for top_k = top_j+1:length(myCommunity)
                 COR(top_j,top_k) = length(intersect(myCommunity{top_j},myCommunity{top_k}))/min(length(myCommunity{top_j}),length(myCommunity{top_k}));
                 if COR(top_j,top_k)>0.4
                     myCommunity{top_j} = union(myCommunity{top_j},myCommunity{top_k});
                     myCommunity{top_k} = [];
                     
                 end
             end
         end
         myCommunity(cellfun(@isempty,myCommunity))=[];
         topTree{treeIndex} = myCommunity;

     maxQ = 0;
     Q = [];
     for tt = 1:length(topTree)
         myLabel = zeros(1,length(AdjMatrix));
         for yy = 1:length(topTree{tt})
            myLabel(topTree{tt}{yy}) = yy;
         end
         Q = [Q,modularity(myLabel,AdjMatrix)];
     end
     [~,sortIndex] = sort(Q,'descend');
     
     [maxQ,maxQ_index] = max(Q);
     
     if length(sortIndex)>=10
         len = 10;
     else
         len = length(sortIndex);
     end
%      myLabel = zeros(len,length(AdjMatrix));

     tic
     tt = 0;
remainNodes2 = [];
     remainNodes3 = [];
    for i = 1:length(remainNodes)
        i_nebor = find(AdjMatrix(remainNodes(i),:)>0);
        i_nebor_label = myLabel(i_nebor);
        i_nebor_label(i_nebor_label==0) =[];
        if ~isempty(i_nebor_label)
            X = unique(i_nebor_label);
            [M,N]=hist(i_nebor_label,X);
            [~,index] = max(M);
            myLabel(remainNodes(i)) = N(index);
        else
%             tt = tt+1;
%            myLabel(remainNodes(i)) = randi(length(myCommunity)); 
        remainNodes2 = [remainNodes2 remainNodes(i)];
        end
    end
    
    for i = 1:length(remainNodes2)
        i_nebor = find(AdjMatrix(remainNodes2(i),:)>0);
        i_nebor_label = myLabel(i_nebor);
        i_nebor_label(i_nebor_label==0) =[];
        if ~isempty(i_nebor_label)
            X = unique(i_nebor_label);
            [M,N]=hist(i_nebor_label,X);
            [~,index] = max(M);
            myLabel(remainNodes2(i)) = N(index);
        else
%             tt = tt+1;
%            myLabel(remainNodes2(i)) = randi(length(myCommunity)); 
         remainNodes3 = [remainNodes3 remainNodes2(i)];
        end
    end
    
    for i = 1:length(remainNodes3)
        i_nebor = find(AdjMatrix(remainNodes3(i),:)>0);
        i_nebor_label = myLabel(i_nebor);
        i_nebor_label(i_nebor_label==0) =[];
        if ~isempty(i_nebor_label)
            X = unique(i_nebor_label);
            [M,N]=hist(i_nebor_label,X);
            [~,index] = max(M);
            myLabel(remainNodes3(i)) = N(index);
        else
            tt = tt+1;
           myLabel(remainNodes3(i)) = randi(length(myCommunity)); 
        end
    end
toc
     

root1 = sprintf('results/%s/remainNodes',name1);
if ~isdir(root1) %判断路径是否存在
    mkdir(root1);
end
 path1=sprintf('results/%s/remainNodes/remainNodes(remainNodesNum=%s)_%s.txt',name1,num2str(length(remainNodes)),num2str(c));
savedata1(path1,remainNodes);
root2 = sprintf('results/%s/addRemainNodesCommunity',name1);
if ~isdir(root2) %判断路径是否存在
    mkdir(root2);
end
 path2=sprintf('results/%s/addRemainNodesCommunity/addRemainNodesCommunity_%s.txt',name1,num2str(c));
savedata1(path2,myLabel);





     modnmi1 = [];
     for i = 1:len
         finalCommunity{i} = topTree{sortIndex(i)};
         qqq = single(modularity(myLabel(i,:),AdjMatrix));%%分派离散点后直接计算Q
         qqq
        modnmi1(i,1) = qqq ;
        if sum(Datalabel) ~= 0
          modnmi1(i,2) = nmi(Datalabel,myLabel(i,:));
       else
          modnmi1(i,2) = 0;
       end
         
         
         
         [myLabel,Q2]=find_error2(single(myLabel),AdjMatrix,qqq);
         
%          allNodes = cell2mat(finalCommunity{i});
%          A=tabulate(allNodes);
%          errorNodes=find(A(:,2)>1);
%          
%           [myLabel,remove,Q2]=find_error(single(myLabel),AdjMatrix,finalCommunity{i},qqq);
%           finalCommunity2 = {};
%           ll = unique(myLabel);
%           for ii = 1:length(ll)
%               finalCommunity2{ii} = find(myLabel==ll(1,ii));
%           end
%            [myLabel,remove,Q2]=find_error(single(myLabel),AdjMatrix,finalCommunity2,Q2);
          
          
          
          
          modnmi = zeros(size(myLabel,1),2);
            for k = 1:size(myLabel,1)
            modnmi(k,1) = (modularity(myLabel(k,:),AdjMatrix)) ;%%纠错后计算
               if sum(Datalabel) ~= 0
                  modnmi(k,2) = nmi(Datalabel,myLabel(k,:));
              else
                  modnmi(k,2) = 0;
              end
           end
            modnmi

     end
     
     clock2 = clock;
     etime(clock2,clock1)
     pathTime=sprintf('results/%s/time/postProcessTime/postProcessTime_%s.txt',name1,num2str(c));
     savedata1(pathTime, etime(clock2,clock1));
    
     
     

   root3 = sprintf('results/%s/findErrorCommunity',name1);
    if ~isdir(root3) %判断路径是否存在
        mkdir(root3);
    end
 path3=sprintf('results/%s/findErrorCommunity/findErrorCommunity_%s.txt',name1,num2str(c));
savedata1(path3,myLabel);
pathResult = sprintf('results/%s/Q_NMI/Q_NMI_%s.txt',name1,num2str(c));
result = [modnmi0;modnmi1;modnmi];
savedata1(pathResult,[result;0 0 ;mean(result(:,1)) mean(result(:,2));max(result(:,1)) max(result(:,2))]);
% modnmi = single(zeros(size(myLabel,1),2));
% for i = 1:size(myLabel,1)
% modnmi(i,1) = single(modularity(myLabel(i,:),AdjMatrix)) ;
% modnmi(i,2) = single(nmi(Datalabel,myLabel(i,:)));
% end
% name

% modnmi
 length(remainNodes)
 tt
% tic
% [myLabel,remove]=find_error(single(myLabel),AdjMatrix,Out_Community,modnmi(:,1));
% 
% modnmi = single(zeros(size(myLabel,1),2));
% for i = 1:size(myLabel,1)
% modnmi(i,1) = single(modularity(myLabel(i,:),AdjMatrix)) ;
% modnmi(i,2) = single(nmi(Datalabel,myLabel(i,:)));
% end
% % name
% modnmi

% qqq
% toc

 end 
 
 
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

function [currentLabel,Q2]=find_error2(currentLabel,AdjMatrix,Q)
    Q2 = Q;
%         M=[];
        
%         numVar=size(AdjMatrix,1);
%             change_node=[];
%             lable= currentLabel;
            %% %%后处理：2种方式-最大的Q lable= ParetoFront1(in,1:numVar);
            
            
            T = 5;
           while T>0
            originalLabel = currentLabel;
            communityNum=unique(currentLabel);
            for i=1:length(communityNum)
                currentCommunity = find(currentLabel==communityNum(1,i));
                if(length(currentCommunity)>=3)
                    for j=1:length(currentCommunity)
%                         m=CLique{i}(j);
%                         k=lable(m);
% k = i;
%                         index_last=find(lable==i);
%                         A=setdiff(index_last,m);
                        neighbors=find(AdjMatrix(currentCommunity(1,j),:));
                        com_max=multi_label(neighbors,currentLabel);
                        currentLabel(1,currentCommunity(1,j))=com_max;
%                         if com_max~=i
%                             change_node=[change_node{in} m];
%                         end
                    end
                end
            end
            QQ = modularity(currentLabel,AdjMatrix);
            if QQ>Q                                                      %有疑问
%                 M=[M in];
               
                Q2 = QQ;
%                 finalLabel = currentLabel;
            else
               currentLabel =  originalLabel;
            end
             T = T-1;
           end
%         error_node=[];
%         for i=1:length(M)
%             error_node=[error_node change_node{M(i)}];
%         end
%         remove=unique(single(error_node));
        
    end

function [ParetoFront1,remove,Q2]=find_error(ParetoFront1,AdjMatrix,CLique,Q)
Q2 = Q;
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
                Q2 = QQ(in);
            end
        end
        error_node=[];
        for i=1:length(M)
            error_node=[error_node change_node{M(i)}];
        end
        remove=unique(single(error_node));
        
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

