%������Ⱥ������Ӧ�ȣ��Բ�ͬ���Ż�Ŀ�꣬�˴���Ҫ��д
%pop_size: ��Ⱥ��С
%chromo_size: Ⱦɫ�峤��

function fitness(pop_size, chromo_size)
global fitness_value;
global pop;
% global G;
global loadCommunity;
global AdjMatrix;
global communityNum_max;
global degree;
global edgeNum;
global modularity
for i=1:pop_size
    fitness_value(i) = 0.;    
end

for i=1:pop_size
    currentCommunity = cell(1,chromo_size);
    for j=1:chromo_size
%         if pop(i,j) == 1
%             fitness_value(i) = fitness_value(i)+2^(j-1);
%         end     
if pop(i,j)==0||length(loadCommunity{j})<pop(i,j)
    pop(i,j) = communityNum_max(1,j);
end
    currentCommunity(1,j) = loadCommunity{j}(1,pop(i,j));
    end
    myLabel = zeros(1,length(AdjMatrix),'single');
    for ii = 1:length(currentCommunity)
        myLabel(currentCommunity{ii}) = ii;
    end
     fitness_value(i) = modularity(myLabel,AdjMatrix,degree,edgeNum);
end

% clear i;
% clear j;
