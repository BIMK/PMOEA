function [a] = karate_test(fid)

        a = 0
    isRealWorld = 1;
    global total_algo_time time_MOEA time_EA
    overlapping = 0;
 Date = {'karate', 'dolphin', 'football', 'netscience_remove', 'CQC', 'blogs', 'hepth', 'hepth1', 'CA-AstroPh_18772_396160', 'CA-CondMat_23133_186936', 'Brightkite_58228', '196591'};
    %Date={'karate','dolphin','football','netscience_remove','CQC','blogs','Hepth','hepth1','CA-AstroPh','CA-CondMat','Brightkite_edges','Gowalla_edges'};
    for i=12:12
        
      save_root = ['result_statistics/',Date{i}];
      if ~isdir(save_root) %判断路径是否存在
           mkdir(save_root);
      end     
      
      path = sprintf('RealWorld/%s.txt',Date{i});
      name=Date{i};
      real_path=sprintf('RealWorldreal_label_%s.txt',Date{i});
      
       t_start = clock;
       
      
        for c=1:1
            
             [AdjMatrix,degree,edgeNum]=pre_process(path,name,real_path,overlapping,isRealWorld,c);
             myLabel = randi([1,2],1,length(AdjMatrix),'single');            
             STR_class = ['AdjMatrix:',class(AdjMatrix),', degree:',class(degree),', edgeNum:', class(edgeNum),', myLabel:', class(myLabel)];
            Q = modularity(myLabel,AdjMatrix,degree,edgeNum);   
            
               myLabel_1 = randi([1,4],1,length(AdjMatrix),'single');               
            Q_1 = modularity(myLabel_1,AdjMatrix,degree,edgeNum); 
            

%             STR = ['Network:',name,', run:',num2str(c),'/20',', this run used time:',num2str(etime(clock,t0)),', total used time:',num2str(etime(clock,t_start))];
            fprintf(fid,'%s\n',STR_class);
            
            
        end
   end
  

end


function  [AdjMatrix,degree,edgeNum]=pre_process(path,name,real_path,overlapping,isRealWorld,c)
    global edgeMatrix 
global AdjMatrix 
global edgeNum
global degree
global total_algo_time time_EA

%%%%%
name1=name;
AdjMatrix = [];
tic;
%% ?????ڽӾ???
%% ?Զ????ݶ?ȡ?????????ݸ?ʽȷ???ڽӾ????????????ǻ???%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hasReal =  exist(real_path,'file');%% ??????ʵ???? ????0????
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
         if AdjMatrix_size(2)>2 %%?ڽӾ?????ʾ
        %     edgeNum = sum(sum(AdjMatrix));
            [edgeMatrix1,edgeMatrix2] = (find(AdjMatrix==1));
            edgeMatrix = [edgeMatrix2,edgeMatrix1];
         else                   %%?߱???ʾ
             AdjMatrix = reIndex(AdjMatrix);        %%?޸?ʱ???ӵģ?????????YHP
             edgeMatrix = AdjMatrix;
             needAddOne = 0;    %%?Ƿ???Ҫ??1
             numVar=(max(max(AdjMatrix(:,1)),max(AdjMatrix(:,2))));
             if find(AdjMatrix==0)>0 %% ??0??ʼ????
                  needAddOne = 1;
                  numVar=numVar+1;
             end
             edgeNum = AdjMatrix_size(1);
              AdjMatrix = Adjreverse(AdjMatrix,numVar,needAddOne);
         end
end

 numVar=single(size(AdjMatrix,1));
 if hasReal >0
     if overlapping == 0 %%???ص?
         Datalabel=(load(real_path));
         if size(Datalabel,2)==2  %%???Ż???Ϊ????--?????š???2????ʽ
             Datalabel=(Datalabel(:,2)');
         end
     else                %%?ص?
         if isRealWorld ==1 %%??ʵ????
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
     if overlapping == 0 %%???ص?
         Datalabel= false(1,numVar);
     else                %%?ص?
         realCommunity = {};
     end
     
     
 end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear EDGE EDGE_Data edgeMatrix1 edgeMatrix2;

degree=single(sum(AdjMatrix,2));
% edgeNum=sum(degree)/2;

%% ?????޸?ʱ???룺Ӧ?Թ?��?????쳣??????Ԥ???���??YHP

index=degree==0;
AdjMatrix(index,:)=[];
AdjMatrix(:,index)=[];
numVar=length(AdjMatrix);

%%  AdjMatrix(logical(eye(size(AdjMatrix)))) = 0;

degree=single(sum(AdjMatrix,2));
edgeNum=sum(degree)/2;


end
