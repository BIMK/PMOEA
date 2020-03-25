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

elitism = true;%ѡ��Ӣ����
pop_size = 100;%��Ⱥ��С
chromo_size = length(community);%Ⱦɫ���С
generation_size = 50;%��������
cross_rate = 0.6;%�������
mutate_rate = 0.01;%�������

[m,~,~,~,bestQ] = GeneticAlgorithm(name,pop_size, chromo_size, generation_size, cross_rate, mutate_rate,elitism,community);
% disp "���Ÿ���"
% m
% disp "������Ӧ��"
% n
% disp "���Ÿ����Ӧ�Ա���ֵ"
% q
% disp "�õ����Ž���Ĵ���"
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