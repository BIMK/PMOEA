%遗传算法主函数
%pop_size: 输入种群大小
%chromo_size: 输入染色体长度
%generation_size: 输入迭代次数
%cross_rate: 输入交叉概率
%cross_rate: 输入变异概率
%elitism: 输入是否精英选择
%m: 输出最佳个体
%n: 输出最佳适应度
%p: 输出最佳个体出现代
%q: 输出最佳个体自变量值

function [m,n,p,q,bestQ] = GeneticAlgorithm(name,pop_size, chromo_size, generation_size, cross_rate, mutate_rate, elitism,community)
tic
global G ; %当前代
global fitness_value;%当前代适应度矩阵
global best_fitness;%历代最佳适应值
global fitness_avg;%历代平均适应值矩阵
global best_individual;%历代最佳个体
global best_generation;%最佳个体出现代
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
initilize(pop_size, chromo_size);%初始化
for G=1:generation_size   
    fitness(pop_size, chromo_size);  %计算适应度 
    rank(pop_size, chromo_size);  %对个体按适应度大小进行排序
    selection(pop_size, chromo_size, elitism);%选择操作
    crossover(pop_size, chromo_size, cross_rate);%交叉操作
    mutation(pop_size, chromo_size, mutate_rate);%变异操作
    best_Q = [best_Q max(fitness_value)];
%     clc;
%     fprintf('%5s网络,已完成%4s%%,耗时%5s秒\n',name,num2str(roundn(G/generation_size*100,-1)),num2str(roundn(toc,-2)));
end
% plotGA(generation_size);%打印算法迭代过程
m = best_individual;%获得最佳个体
n = best_fitness;%获得最佳适应度
p = best_generation;%获得最佳个体出现代
bestQ = max(best_Q);
%获得最佳个体变量值，对不同的优化目标，此处需要改写
q = 0.;
% for j=1:chromo_size
%     if best_individual(j) == 1
%             q = q+2^(j-1);
%     end 
% end
% q = -1+q*(3.-(-1.))/(2^chromo_size-1);

% clear i;
% clear j;