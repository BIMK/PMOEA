function MAIN_START(path,name,real_path,max_gen)

% clc;
% clear;
tic;

global numObjectives popsize numVar  max_gen mutate_posibility M niche adj_mat V ll pc pg pm
global   strNetwork Datalable idealp weights neighbors  edgeslist  degree
global Thirty_Run_maxQ Thirty_Run_maxNMI   Adj_mat Node AdjMatrix Node3 matrix Adj_mat

%max_gen         =400;
pg=1;

M = 2;
pm=0.3;
pc=1;
runtimes=1;
%% 多实验
numVar=5000;
EDGE=load(path);
Datalable=load(real_path);
Datalable=Datalable(:,2)';
AdjMatrix=Adjreverse(EDGE,numVar)
strNetwork=name;
%AdjMatrix = load(path);
%Datalable=load(real_path);
%Datalable=Datalable';
%% 公茂果团队的真实网络和GNEtend网络
% strNetwork='0.45';
% AdjMatrix=load('GNExtend/0.45.txt');
% Datalable=load('GNExtend/real0.45.txt');
% AdjMatrix = load('RealWorld\football.txt');
% Datalable=load('RealWorld\real_label_football.txt');
%% 人工随机网络
% true_partion([1:32])=1; true_partion([33:64])=2;true_partion([65:96])=3;true_partion([97:128])=4;%%
% Datalable=true_partion;
% DATAfile = '.\基准随机网络\Pin0.4_0.txt';
% AdjMatrix=vary_edge_Adj(DATAfile);
%% 获取合并后的团及相应的信息
[Node1,matrix,degree1,edgeslist]=Get_Cliques(AdjMatrix,Datalable);

%% 文件数据保存及进化得到parate面
root = sprintf('results1/%s/ParetoFront',strNetwork);  %%建个文件保存实验数据
if ~isdir(root) %判断路径是否存在
    mkdir(root); 
end
root = sprintf('results1/%s/metrics',strNetwork);
if ~isdir(root) %判断路径是否存在
    mkdir(root);
end
Thirty_Run_maxQ=[];
Thirty_Run_maxNMI=[];
for mg=1:runtimes
    fprintf('第%s次运行\n',num2str(mg));
    
    TMOEAD(matrix,Node1,mg,degree1);
    fprintf('\n');
  
end 
[~,index1]=max(Thirty_Run_maxQ(:,1));%找到最大的Q索引
[~,index2]=max(Thirty_Run_maxNMI(:,2));%找到最大的NMI索引
fprintf('Thirtyrun_results:\n');
fprintf('maxQ =    %g %g\n',Thirty_Run_maxQ(index1,1),Thirty_Run_maxQ(index1,2));    %输出最大的Q时的Q和NMI
fprintf('maxNMI =  %g %g\n',Thirty_Run_maxNMI(index2,1),Thirty_Run_maxNMI(index2,2));%输出最大的NMI时的NMI和Q
path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxQ.txt',strNetwork,strNetwork);%保存实验数据取平均值
savedata1(path,[Thirty_Run_maxQ;0 0;mean(Thirty_Run_maxQ(:,1)) mean(Thirty_Run_maxQ(:,2))]); 
path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxNMI.txt',strNetwork,strNetwork);
savedata1(path,[Thirty_Run_maxNMI;0 0;mean(Thirty_Run_maxNMI(:,1)) mean(Thirty_Run_maxNMI(:,2))]);
clc;

clear;
end