%单点交叉操作
%pop_size: 种群大小
%chromo_size: 染色体长度
%cross_rate: 交叉概率

function crossover(pop_size, chromo_size, cross_rate)
global pop;
global communityNum_max;
for i=1:2:pop_size
    if(rand < cross_rate)
        cross_pos = round(rand * chromo_size);
        if or (cross_pos == 0, cross_pos == 1)
            continue;
        end
        for j=cross_pos:chromo_size
            temp = pop(i,j);
            if pop(i,j) >= pop(i+1,j)%%防止越界
                pop(i,j) = pop(i+1,j);
                pop(i+1,j) = mod(temp,communityNum_max(1,j));
%                 if pop(i+1,j) == 0
%                    pop(i+1,j) =  pop(i,communityNum_max(1,j)); 
% %                    pop(i,j) = temp;
%                    
%                 end
pop(i,j) = temp;
            else
                pop(i,j) = mod(pop(i+1,j),communityNum_max(1,j));
%                 if pop(i,j) == 0
%                    pop(i,j) =  pop(i,communityNum_max(1,j)); 
%                    pop(i+1,j) = temp;
%                    
%                 end
pop(i+1,j) = temp;
            end
            
            if pop(i,j)==0
                pop(i,j) = communityNum_max(1,j);
            end
            
            if pop(i+1,j)==0
                pop(i+1,j) = communityNum_max(1,j);
            end
            
        end
    end
end


% clear i;
% clear j;
% clear temp;
% clear cross_pos;