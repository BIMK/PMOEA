function  PMOEA_EA(path,name,real_path,overlapping,isRealWorld,c)
%%%%%
global edgeMatrix 
global AdjMatrix 
global edgeNum
global degree
global total_algo_time time_EA
global modularity

%%%%%
name1=name;
AdjMatrix = [];
tic;
%% �����ڽӾ���
%% �Զ����ݶ�ȡ���������ݸ�ʽȷ���ڽӾ�����������ǻ���%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hasReal =  exist(real_path,'file');%% ������ʵ���� ����0����
 networkData = sprintf('adjMatrix_coreNodes/%s/%s.mat',name1,name1);
 hasNetWorkAdj =  exist(networkData,'file');
if hasNetWorkAdj
    load(networkData);
else
    AdjMatrix = single(load(path));
end
networkData2 = sprintf('adjMatrix_coreNodes/%s/coreNodes.mat',name1);
load(networkData2);
if isempty(AdjMatrix)||size(AdjMatrix,2) == 2
%      AdjMatrix = single(load(path));
         AdjMatrix_size = size(AdjMatrix);
         if AdjMatrix_size(2)>2 %%�ڽӾ����ʾ
        %     edgeNum = sum(sum(AdjMatrix));
            [edgeMatrix1,edgeMatrix2] = (find(AdjMatrix==1));
            edgeMatrix = [edgeMatrix2,edgeMatrix1];
         else                   %%�߱��ʾ
             AdjMatrix = reIndex(AdjMatrix);        %%�޸�ʱ��ӵģ��������YHP
             edgeMatrix = AdjMatrix;
             needAddOne = 0;    %%�Ƿ���Ҫ��1
             numVar=(max(max(AdjMatrix(:,1)),max(AdjMatrix(:,2))));
             if find(AdjMatrix==0)>0 %% ��0��ʼ���
                  needAddOne = 1;
                  numVar=numVar+1;
             end
             edgeNum = AdjMatrix_size(1);
              AdjMatrix = Adjreverse(AdjMatrix,numVar,needAddOne);
         end
end

if ~isa(AdjMatrix,'logical')
    AdjMatrix = logical(AdjMatrix);
end

 numVar=single(size(AdjMatrix,1));
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
 clear EDGE EDGE_Data edgeMatrix1 edgeMatrix2;

degree=single(sum(AdjMatrix,2));
% edgeNum=sum(degree)/2;

%% ����޸�ʱ���룺Ӧ�Թ�������쳣�����Ԥ������YHP

index=degree==0;
AdjMatrix(index,:)=[];
AdjMatrix(:,index)=[];
numVar=length(AdjMatrix);

%%  AdjMatrix(logical(eye(size(AdjMatrix)))) = 0;

degree=single(sum(AdjMatrix,2));
edgeNum=sum(degree)/2;

%% ��ȡ�� MOEA ����Ľ⼯��Ӧ������MOEA������⼯����ͬһ��Ŀ¼�£�
community = {};
networkDataRoot = sprintf('results/%s/coreCommunityFirstFront/coreCommunity_%s/',name1,num2str(c));%%��ȡ��c��MOEA���е����ż���
    path = networkDataRoot;  
    fileExt = '*.mat';  
    files = dir(fullfile(path,fileExt));  
    for i=1:size(files,1)  
        fileName = strcat(path,files(i,1).name);  
        coreCommunity = load(fileName);
%         community = [community {coreCommunity.chromosomes}];
community = coreCommunity.cellCommunity;
    end; 
    
     pathTime = sprintf('results/%s/time/postProcessTime/',name1);
    if ~isdir(pathTime) %�ж�·���Ƿ����
        mkdir(pathTime);
    end
    pathResult = sprintf('results/%s/Q_NMI/',name1);
    if ~isdir(pathResult) %�ж�·���Ƿ����
        mkdir(pathResult);
    end
     pathResultError = sprintf('results/%s/errorNodes/',name1);
    if ~isdir(pathResultError) %�ж�·���Ƿ����
        mkdir(pathResultError);
    end
    
      pathResultByEA = sprintf('results/%s/CommunityByEA/',name1);
    if ~isdir(pathResultByEA) %�ж�·���Ƿ����
        mkdir(pathResultByEA);
    end
    
    clock1 = clock;
    %% ������MOEA����ĸ���⼯��ѡ��
    [communityIndex,bestQ] = run_EA(name,community);

    findNodes = [];
    for ii = 1:length(communityIndex)
        myCommunity{ii} = community{ii}{communityIndex(ii)}(1,1:end-2);
        findNodes = unique([findNodes myCommunity{ii}]);
    end
    save([pathResultByEA 'CommunityByEA.mat'],'myCommunity', '-v7.3');
 
    myLabel = zeros(1,length(Datalabel),'single');
    for ii = 1:length(myCommunity)
        myLabel(myCommunity{ii}) = ii;
    end
    modnmi0 = [];
    modnmi0(1,1) = bestQ ;%%��������ɢ��ֱ�Ӽ���
    if sum(Datalabel) ~= 0
      modnmi0(1,2) = nmi(Datalabel,myLabel(1,:));
    else
      modnmi0(1,2) = 0;
    end

remainNodes = setdiff(1:length(AdjMatrix),findNodes);

%%%%%%%%%%%%%%%%���ƶȴ�������ں�%%%%%%%%%%%
topTree = {};
 finalCommunity = {};
     treeIndex = 1;
     topTree{treeIndex} = myCommunity;

         COR = zeros(length(myCommunity),length(myCommunity),'single');
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
         myLabel = zeros(1,length(AdjMatrix),'single');
         for yy = 1:length(topTree{tt})
            myLabel(topTree{tt}{yy}) = yy;
         end
         Q = [Q,modularity(myLabel,AdjMatrix,degree,edgeNum)];
     end
     [~,sortIndex] = sort(Q,'descend');
     
     [maxQ,maxQ_index] = max(Q);
     
     if length(sortIndex)>=10
         len = 10;
     else
         len = length(sortIndex);
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
%%%%%%%%%%%%%%%%ʣ��ڵ㴦��%%%%%%%%
%      tic
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
    % toc
    %%%%%%%%%%%%%%%%ʣ��ڵ㴦��%%%%%%%%%%%

     time_EA = toc
     total_algo_time = total_algo_time + time_EA;

root1 = sprintf('results/%s/remainNodes',name1);
if ~isdir(root1) %�ж�·���Ƿ����
    mkdir(root1);
end
 path1=sprintf('results/%s/remainNodes/remainNodes(remainNodesNum=%s)_%s.txt',name1,num2str(length(remainNodes)),num2str(c));
savedata1(path1,remainNodes);
root2 = sprintf('results/%s/addRemainNodesCommunity',name1);
if ~isdir(root2) %�ж�·���Ƿ����
    mkdir(root2);
end
 path2=sprintf('results/%s/addRemainNodesCommunity/addRemainNodesCommunity_%s.txt',name1,num2str(c));
savedata1(path2,myLabel);

     modnmi1 = [];
     for i = 1:len
         finalCommunity{i} = topTree{sortIndex(i)};
         Q0 = single(modularity(myLabel(i,:),AdjMatrix,degree,edgeNum));%%������ɢ���
        modnmi1(i,1) = Q0 ;
        if sum(Datalabel) ~= 0
          modnmi1(i,2) = nmi(Datalabel,myLabel(i,:));
       else
          modnmi1(i,2) = 0;
       end
         
% % % %          allNodes = cell2mat(finalCommunity{i});
% % % %          A=tabulate(allNodes);
% % % %          errorNodes=find(A(:,2)>1);
         %%����
           [myLabel,Q2]=find_error2(single(myLabel),AdjMatrix,Q0,degree,edgeNum);
          
          modnmi = zeros(size(myLabel,1),2,'single');
            for k = 1:size(myLabel,1)
            modnmi(k,1) = (modularity(myLabel(k,:),AdjMatrix,degree,edgeNum)) ;%%
               if sum(Datalabel) ~= 0
                  modnmi(k,2) = nmi(Datalabel,myLabel(k,:));
              else
                  modnmi(k,2) = 0;
              end
           end
            modnmi

     end
     
     clock2 = clock;
     etime(clock2,clock1);
     pathTime=sprintf('results/%s/time/postProcessTime/postProcessTime_%s.txt',name1,num2str(c));
     savedata1(pathTime, etime(clock2,clock1));
     
    
% % % %      labelIndex = unique(myLabel);
% % % %      finalCommunityNum = length(labelIndex);
% % % %      COMMUNITY = {};
% % % %      for iii = 1:finalCommunityNum
% % % %          COMMUNITY{iii} = find(myLabel == labelIndex(iii));
% % % %      end
     
% % % %       afterErrorCorrectionErrorNodes = zeros(1,finalCommunityNum,'single');
% % % %     for oo = 1:finalCommunityNum
% % % %         afterErrorCorrectionErrorNodes(oo) = length(unique(Datalabel(1,COMMUNITY{oo})))-1;%%����ʵ������������
% % % %     end
% % % %     
% % % %      avgErrorNodes = sum(afterErrorCorrectionErrorNodes)/finalCommunityNum;%%ÿ���ŵĴ���ڵ���
% % % %      
% % % %      %%���д���ڵ���������������
% % % %      pathResultError2=sprintf('results/%s/errorNodes/afterErrorCorrectionErrorNodes(%s_%s=%s)_%s.txt',name1,num2str(sum(afterErrorCorrectionErrorNodes)),num2str(finalCommunityNum),num2str(avgErrorNodes),num2str(c));
% % % %      savedata1(pathResultError2,afterErrorCorrectionErrorNodes);
     
     
   root3 = sprintf('results/%s/findErrorCommunity',name1);
    if ~isdir(root3) %�ж�·���Ƿ����
        mkdir(root3);
    end
 path3=sprintf('results/%s/findErrorCommunity/findErrorCommunity(cNum=%s)_%s.txt',name1,num2str(length(unique(myLabel))),num2str(c));
savedata1(path3,myLabel);
pathResult = sprintf('results/%s/Q_NMI/Q_NMI_%s.txt',name1,num2str(c));
result = [modnmi0;modnmi1;modnmi];
savedata1(pathResult,[result;0 0 ;mean(result(:,1)) mean(result(:,2));max(result(:,1)) max(result(:,2))]);

save_root = ['result_statistics/',name,'/Q_NMI'];
    if ~isdir(save_root) %�ж�·���Ƿ����
             mkdir(save_root);    
    end
dlmwrite([save_root,'/Q_NMI_',num2str(c),'.txt'], [result;0 0 ;mean(result(:,1)) mean(result(:,2));max(result(:,1)) max(result(:,2))],' ')

 end 
 


function [currentLabel,Q2]=find_error2(currentLabel,AdjMatrix,Q,degree,edgeNum)
    Q2 = Q;

            %% %%����2�ַ�ʽ-����Q lable= ParetoFront1(in,1:numVar);
            
            
            T = 5;
           while T>0
            originalLabel = currentLabel;
            communityNum=unique(currentLabel);
            for i=1:length(communityNum)
                currentCommunity = find(currentLabel==communityNum(1,i));
                if(length(currentCommunity)>=3)
                    for j=1:length(currentCommunity)

                        neighbors=find(AdjMatrix(currentCommunity(1,j),:));
                        if ~isempty(neighbors)      %% ���ݱ���ӵ�if����YHP
                            com_max=multi_label(neighbors,currentLabel);
                            currentLabel(1,currentCommunity(1,j))=com_max;
                        end
                    end
                end
            end
            QQ = modularity(currentLabel,AdjMatrix,degree,edgeNum);
            if QQ>Q                                                   
               
                Q2 = QQ;
            else
               currentLabel =  originalLabel;
            end
             T = T-1;
           end

        
    end



