function last_operate(string_m,AdjacentMatrix,Datalable,e,degree1,Node3)

global  Thirty_Run_maxQ Thirty_Run_maxNMI

strNetwork=string_m;
path = sprintf('results1/%s/ParetoFront/MODPSO1_%s_PF%d.txt',string_m,string_m,e);
ParetoFront=load(path);

% AdjacentMatrix = load('GNExtend\0.45.txt');
% Datalable=load('GNExtend\real0.45.txt');
% AdjacentMatrix = load('RealWorld\football.txt');
% Datalable=load('RealWorld\real_label_football.txt');
Adj_mat=AdjacentMatrix;
% m=find(B==1);
% n=find(B==4);
% B(m)=B(m)+3;
% B(n)=B(n)-3;
% V=Datalable-B;
% index=find(V~=0);
%[Node1,matrix,degree1,edgeslist]=Get_Cliques(Adj_mat,Datalable);
numVar=size(Adj_mat,1);

V=length(degree1);

for in=1:size(ParetoFront,1)
    A=ParetoFront(in,1:V);
   ParetoFront1(in,1:numVar)= decode1(A,Node3);
   change_node{in}=[];
    %% %%后处理：2种方式-最大的Q
    lable= ParetoFront1(in,1:numVar);
    for i=1:length(degree1)
      if(length(Node3{i})>=1)
       for j=1:length(Node3{i})
           
           m=Node3{i}(j);
         if m==114
             wo=1;
         end
           k=lable(m);
           index_last=find(lable==k);
           A=setdiff(index_last,m);
%            if(fff(A,Adj_mat,ll)>fff(index_last,Adj_mat,ll))
               neighbors=find(Adj_mat(m,:));
            com_max=Max(neighbors,lable);
               ParetoFront1(in,m)=com_max;
               if com_max~=k
                   change_node{in}=[change_node{in} m];
               end
       
          
       end
      end
    end
           
         
end
for i=1:size(ParetoFront1,1)
    Q(i,1) = Modularity(ParetoFront1(i,1:numVar),AdjacentMatrix); %%计算模块度
    NMI(i,1) = nmi(ParetoFront1(i,1:numVar),Datalable);%%计算精度
end


metrics=[Q NMI];
[~,index1]=max(Q);
[~,index2]=max(NMI);
fprintf('Paretosize = %g\n',size(ParetoFront1,1));
fprintf('maxQ :    %g  %g\n',metrics(index1,1),metrics(index1,2));
fprintf('maxNMI :  %g  %g\n',metrics(index2,1),metrics(index2,2));
path = sprintf('results1/%s/metrics/MODPSO1_%s_metrics%d.txt',strNetwork,strNetwork,e);
savedata1(path,metrics);
Thirty_Run_maxQ=[Thirty_Run_maxQ;metrics(index1,:)];
Thirty_Run_maxNMI=[Thirty_Run_maxNMI;metrics(index2,:)];
 end
