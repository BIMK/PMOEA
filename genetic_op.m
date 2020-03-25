function [ f ] = genetic_op( parent_chromosome,chromosomes,AdjMatrix,pc,pm,V,M,ll )
%GENETIC_OP Summary of this function goes here
[N,m] = size(parent_chromosome);
[N1,m1]=size(chromosomes);
r0=0;
N=10;
clear m m1
% crossover
p = 1;

for i = 1 : round(N/2) 
    % Select the first parent
    parent_1 = ceil(N*rand(1));  
    r0=rand(1);
  if r0<=0.9
    % Select the second parent
    parent_2 = ceil(N*rand(1));     
    % Make sure both the parents are not the same. 
    while parent_2 == parent_1
        parent_2 = ceil(N*rand(1));
    end
       parent_1 = parent_chromosome(parent_1,:);
       parent_2 = parent_chromosome(parent_2,:);
  else
       parent_2=N+floor((N1-N)*rand(1));
       parent_1 = parent_chromosome(parent_1,:);
       parent_2 = chromosomes(parent_2,:);
  end
    % Get the chromosome information for each randomnly selected parents
   % parent_1 = parent_chromosome(parent_1,:);
   % parent_2 = parent_chromosome(parent_2,:);
    if rand<=pc
        % Perform corssover for each decision variable in the chromosome.
        child_1 = parent_1(1:V);
        child_2 = parent_2(1:V);
        position1 = ceil(V*rand(1));
        position2 = ceil(V*rand(1));
%         while position2 == position1
%             position2 = ceil(V*rand(1));
%         end
        minp = min(position1,position2);
        maxp = max(position1,position2);
       r1=0;
        r1=floor(3*rand(1)+1);   %%floor（A）取A最小整数
        %%分三段随机交叉
        if r1==1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = child_1(1:minp);
           child_1(1:minp) = child_2(1:minp);
            child_2(1:minp) = temp;
        % Evaluate the objective function for the offsprings and as before
        % concatenate the offspring chromosome with objective value.
        child_1(:,V + 1: M + V) = evaluate_objective(child_1,AdjMatrix,ll);
        child_2(:,V + 1: M + V) = evaluate_objective(child_2,AdjMatrix,ll);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       elseif r1==2
            temp = child_1(minp:maxp);
            child_1(minp:maxp) = child_2(minp:maxp);
            child_2(minp:maxp) = temp;
            % Evaluate the objective function for the offsprings and as before
            % concatenate the offspring chromosome with objective value.
            child_1(:,V + 1: M + V) = evaluate_objective(child_1,AdjMatrix,ll);
            child_2(:,V + 1: M + V) = evaluate_objective(child_2,AdjMatrix,ll);  
        elseif r1==3
            temp = child_1(maxp:V);
            child_1(maxp:V) = child_2(maxp:V);
            child_2(maxp:V) = temp;
            % Evaluate the objective function for the offsprings and as before
            % concatenate the offspring chromosome with objective value.
           child_1(:,V + 1: M + V) = evaluate_objective(child_1,AdjMatrix,ll);
            child_2(:,V + 1: M + V) = evaluate_objective(child_2,AdjMatrix,ll);  
        end
    else
        child_1 = parent_1(1:V+M);
        child_2 = parent_2(1:V+M);  
    end
    f(p,:) = child_1;
    f(p+1,:) = child_2;
    p = p + 2;
end

% mutation
degree=sum(AdjMatrix,1);
[NN,m] = size(f);
clear m
for i = 1 : NN
    for j = 1:V
        if rand<=pm
            if degree(j) >= 2  %%邻接点有2个以上
                temp = f(i,j);
                while 1
                    neighbours=find(AdjMatrix(j,:)==1);
                    f(i,j) = neighbours(ceil(rand*degree(j)));
                    if temp ~= f(i,j)
                        break;
                    end
                end         
            end     
        end
    end
    f(i,V + 1: M + V) = evaluate_objective(f(i,1:V),AdjMatrix,ll);
end
end

