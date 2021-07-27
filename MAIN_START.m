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
%% ��ʵ��
numVar=5000;
EDGE=load(path);
Datalable=load(real_path);
Datalable=Datalable(:,2)';
AdjMatrix=Adjreverse(EDGE,numVar)
strNetwork=name;
%AdjMatrix = load(path);
%Datalable=load(real_path);
%Datalable=Datalable';
%% ��ï���Ŷӵ���ʵ�����GNEtend����
% strNetwork='0.45';
% AdjMatrix=load('GNExtend/0.45.txt');
% Datalable=load('GNExtend/real0.45.txt');
% AdjMatrix = load('RealWorld\football.txt');
% Datalable=load('RealWorld\real_label_football.txt');
%% �˹��������
% true_partion([1:32])=1; true_partion([33:64])=2;true_partion([65:96])=3;true_partion([97:128])=4;%%
% Datalable=true_partion;
% DATAfile = '.\��׼�������\Pin0.4_0.txt';
% AdjMatrix=vary_edge_Adj(DATAfile);
%% ��ȡ�ϲ�����ż���Ӧ����Ϣ
[Node1,matrix,degree1,edgeslist]=Get_Cliques(AdjMatrix,Datalable);

%% �ļ����ݱ��漰�����õ�parate��
root = sprintf('results1/%s/ParetoFront',strNetwork);  %%�����ļ�����ʵ������
if ~isdir(root) %�ж�·���Ƿ����
    mkdir(root); 
end
root = sprintf('results1/%s/metrics',strNetwork);
if ~isdir(root) %�ж�·���Ƿ����
    mkdir(root);
end
Thirty_Run_maxQ=[];
Thirty_Run_maxNMI=[];
for mg=1:runtimes
    fprintf('��%s������\n',num2str(mg));
    
    TMOEAD(matrix,Node1,mg,degree1);
    fprintf('\n');
  
end 
[~,index1]=max(Thirty_Run_maxQ(:,1));%�ҵ�����Q����
[~,index2]=max(Thirty_Run_maxNMI(:,2));%�ҵ�����NMI����
fprintf('Thirtyrun_results:\n');
fprintf('maxQ =    %g %g\n',Thirty_Run_maxQ(index1,1),Thirty_Run_maxQ(index1,2));    %�������Qʱ��Q��NMI
fprintf('maxNMI =  %g %g\n',Thirty_Run_maxNMI(index2,1),Thirty_Run_maxNMI(index2,2));%�������NMIʱ��NMI��Q
path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxQ.txt',strNetwork,strNetwork);%����ʵ������ȡƽ��ֵ
savedata1(path,[Thirty_Run_maxQ;0 0;mean(Thirty_Run_maxQ(:,1)) mean(Thirty_Run_maxQ(:,2))]); 
path = sprintf('results1/%s/metrics/MODPSO1_%s_Thirty_Run_maxNMI.txt',strNetwork,strNetwork);
savedata1(path,[Thirty_Run_maxNMI;0 0;mean(Thirty_Run_maxNMI(:,1)) mean(Thirty_Run_maxNMI(:,2))]);
clc;

clear;
end