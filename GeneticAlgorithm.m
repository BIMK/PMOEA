%�Ŵ��㷨������
%pop_size: ������Ⱥ��С
%chromo_size: ����Ⱦɫ�峤��
%generation_size: �����������
%cross_rate: ���뽻�����
%cross_rate: ����������
%elitism: �����Ƿ�Ӣѡ��
%m: �����Ѹ���
%n: ��������Ӧ��
%p: �����Ѹ�����ִ�
%q: �����Ѹ����Ա���ֵ

function [m,n,p,q,bestQ] = GeneticAlgorithm(name,pop_size, chromo_size, generation_size, cross_rate, mutate_rate, elitism,community)
tic
global G ; %��ǰ��
global fitness_value;%��ǰ����Ӧ�Ⱦ���
global best_fitness;%���������Ӧֵ
global fitness_avg;%����ƽ����Ӧֵ����
global best_individual;%������Ѹ���
global best_generation;%��Ѹ�����ִ�
global loadCommunity;
% global AdjMatrix;
global best_Q;
% global degree;
% global edgeNum;
% degree=single(sum(AdjMatrix,2));
% edgeNum=sum(degree)/2;

% best_individual = 
best_Q =[];
% networkData = sprintf('c:/adjMatrix_coreNodes/%s/%s.mat',name,name);
% load(networkData);
loadCommunity = community;

fitness_avg = zeros(generation_size,1,'single');

fitness_value(pop_size) = 0.;
best_fitness = 0.;
best_generation = 0;
initilize(pop_size, chromo_size);%��ʼ��
for G=1:generation_size   
    fitness(pop_size, chromo_size);  %������Ӧ�� 
    rank(pop_size, chromo_size);  %�Ը��尴��Ӧ�ȴ�С��������
    selection(pop_size, chromo_size, elitism);%ѡ�����
    crossover(pop_size, chromo_size, cross_rate);%�������
    mutation(pop_size, chromo_size, mutate_rate);%�������
    best_Q = [best_Q max(fitness_value)];
%     clc;
%     fprintf('%5s����,�����%4s%%,��ʱ%5s��\n',name,num2str(roundn(G/generation_size*100,-1)),num2str(roundn(toc,-2)));
end
% plotGA(generation_size);%��ӡ�㷨��������
m = best_individual;%�����Ѹ���
n = best_fitness;%��������Ӧ��
p = best_generation;%�����Ѹ�����ִ�
bestQ = max(best_Q);
%�����Ѹ������ֵ���Բ�ͬ���Ż�Ŀ�꣬�˴���Ҫ��д
q = 0.;
% for j=1:chromo_size
%     if best_individual(j) == 1
%             q = q+2^(j-1);
%     end 
% end
% q = -1+q*(3.-(-1.))/(2^chromo_size-1);

% clear i;
% clear j;