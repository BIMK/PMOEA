%����������
%pop_size: ��Ⱥ��С
%chromo_size: Ⱦɫ�峤��
%cross_rate: �������
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