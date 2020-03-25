function [communityIndex,bestQ] = run_ga(name,community)
% global name;
global communityNum_max;
% global AdjMatrix;
% name = 'CQC';
% communityNum_max = [];
%     networkDataRoot = sprintf('c:/adjMatrix_coreNodes/%s/coreCommunity/',name);
%     path = networkDataRoot;  
%     fileExt = '*.mat';  
%     files = dir(fullfile(path,fileExt));  
%     for i=1:size(files,1)  
%         fileName = strcat(path,files(i,1).name);  
%         coreCommunity = load(fileName);
%         community = [community {coreCommunity.chromosomes}];
% %         nodesLength = [];
% %         for jj = 1:length(coreCommunity.chromosomes)
% %             nodesLength(jj)  = length(coreCommunity.chromosomes{jj})-2;
% %         end
% %         [~,maxLenIndex] = max(nodesLength);
% %         communityNum_max = [communityNum_max [length(coreCommunity.chromosomes);maxLenIndex]];
%     end; 
    
community = cellfun(@(x)delete_obj(x),community,'UniformOutput',false);
communityNum_max = cellfun('length',community);

elitism = true;%选择精英操作
pop_size = 100;%种群大小
chromo_size = length(community);%染色体大小
generation_size = 50;%迭代次数
cross_rate = 0.6;%交叉概率
mutate_rate = 0.01;%变异概率

[m,~,~,~,bestQ] = GeneticAlgorithm(name,pop_size, chromo_size, generation_size, cross_rate, mutate_rate,elitism,community);
% disp "最优个体"
% m
% disp "最优适应度"
% n
% disp "最优个体对应自变量值"
% q
% disp "得到最优结果的代数"
% p
communityIndex = m;


end
function cell = delete_obj(cell)
for ci = 1:length(cell)
    if length(cell{ci})>2

    cell{ci}(end-1:end) =[];
    end
end
end