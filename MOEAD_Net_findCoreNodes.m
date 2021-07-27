function  MOEAD_Net_findCoreNodes(path,name,real_path,overlapping,isRealWorld)
%%%%%
global edgeMatrix 
global AdjMatrix 
global coreNodes
global degree
global avgDegree


name1=name;
tic
%% 导入邻接矩阵
%% 自动根据读取的网络数据格式确定邻接矩阵和有无真是划分%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hasReal =  exist(real_path,'file');%% 有无真实划分 大于0，有
 AdjMatrix = (load(path));
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
      AdjMatrix=(Adjreverse(AdjMatrix,numVar,needAddOne));
 end


%  numVar=single(size(AdjMatrix,1));
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
AdjMatrix = single(AdjMatrix);
 clear EDGE EDGE_Data edgeMatrix1 edgeMatrix2;
 
 degree=sum(AdjMatrix,1);
 index=find(degree==0);
 AdjMatrix(index,:)=[];
 AdjMatrix(:,index)=[];
 
degree = single(sum(AdjMatrix));
avgDegree = floor(sum(degree)/length(AdjMatrix));

 
strNetwork=name;
networkDataRoot = sprintf('adjMatrix_coreNodes/%s/',name1);
if ~isdir(networkDataRoot) %判断路径是否存在
    mkdir(networkDataRoot);
end

save([networkDataRoot, name1], 'AdjMatrix');

coreNodes = single(sort(findCoreNodesByDegree(AdjMatrix))); 
coreNodes
save([networkDataRoot, 'coreNodes'], 'coreNodes');

end