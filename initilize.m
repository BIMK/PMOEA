%��ʼ����Ⱥ
%pop_size: ��Ⱥ��С
%chromo_size: Ⱦɫ�峤��

function initilize(pop_size, chromo_size)
global pop;
global loadCommunity;
pop = zeros(pop_size,chromo_size);
global best_individual;%������Ѹ���
best_individual = zeros(1,chromo_size);
for i=1:pop_size
    for j=1:chromo_size
        pop(i,j) = randi(length(loadCommunity{j}));
    end
end
clear i;
clear j;