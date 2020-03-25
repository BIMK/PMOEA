%初始化种群
%pop_size: 种群大小
%chromo_size: 染色体长度

function initilize(pop_size, chromo_size)
global pop;
global loadCommunity;
pop = zeros(pop_size,chromo_size);
global best_individual;%历代最佳个体
best_individual = zeros(1,chromo_size);
for i=1:pop_size
    for j=1:chromo_size
        pop(i,j) = randi(length(loadCommunity{j}));
    end
end
clear i;
clear j;