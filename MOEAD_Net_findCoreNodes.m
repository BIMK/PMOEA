function  MOEAD_Net_findCoreNodes(path,name,real_path,overlapping,isRealWorld)
%%%%%
global edgeMatrix 
global AdjMatrix 
global coreNodes
global degree
global avgDegree


name1=name;
tic
%% �����ڽӾ���
%% �Զ����ݶ�ȡ���������ݸ�ʽȷ���ڽӾ�����������ǻ���%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hasReal =  exist(real_path,'file');%% ������ʵ���� ����0����
 AdjMatrix = (load(path));
 AdjMatrix_size = size(AdjMatrix);
 if AdjMatrix_size(2)>2 %%�ڽӾ����ʾ
%     edgeNum = sum(sum(AdjMatrix));
    [edgeMatrix1,edgeMatrix2] = (find(AdjMatrix==1));
    edgeMatrix = [edgeMatrix2,edgeMatrix1];
 else                   %%�߱��ʾ
     edgeMatrix = AdjMatrix;
     needAddOne = 0;    %%�Ƿ���Ҫ��1
     numVar=(max(max(AdjMatrix(:,1)),max(AdjMatrix(:,2))));
     if find(AdjMatrix==0)>0 %% ��0��ʼ���
          needAddOne = 1;
          numVar=numVar+1;
     end
     edgeNum = AdjMatrix_size(1);
      AdjMatrix=(Adjreverse(AdjMatrix,numVar,needAddOne));
 end


%  numVar=single(size(AdjMatrix,1));
 if hasReal >0
     if overlapping == 0 %%���ص�
         Datalabel=(load(real_path));
         if size(Datalabel,2)==2  %%���Ż���Ϊ����--�����š���2����ʽ
             Datalabel=(Datalabel(:,2)');
         end
     else                %%�ص�
         if isRealWorld ==1 %%��ʵ����
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
     if overlapping == 0 %%���ص�
         Datalabel= false(1,numVar);
     else                %%�ص�
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
if ~isdir(networkDataRoot) %�ж�·���Ƿ����
    mkdir(networkDataRoot);
end

save([networkDataRoot, name1], 'AdjMatrix');

coreNodes = single(sort(findCoreNodesByDegree(AdjMatrix))); 
coreNodes
save([networkDataRoot, 'coreNodes'], 'coreNodes');

end