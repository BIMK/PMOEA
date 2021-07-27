%单点变异操作
%pop_size: 种群大小
%chromo_size: 染色体长度
%cross_rate: 变异概率
function mutation(pop_size, chromo_size, mutate_rate)
global pop;
global communityNum_max;
for i=1:pop_size
    if rand < mutate_rate
%         mutate_pos = round(rand*chromo_size);
        mutate_pos = randi(chromo_size);
%         if mutate_pos == 0
%             continue;
%         end
        if communityNum_max(1,mutate_pos)>1
            randomIndex = 1:communityNum_max(1,mutate_pos);
            randomIndex = setdiff(randomIndex,pop(i, mutate_pos));
            randomIndex(randi(length(randomIndex)));
            pop(i,mutate_pos);
             randomIndex(randi(length(randomIndex)));
%             pop(i,mutate_pos) = 1 - pop(i, mutate_pos);
            pop(i,mutate_pos) =  randomIndex(randi(length(randomIndex)));
        else
            continue;
        end
    end
end

% clear i;
% clear mutate_pos;