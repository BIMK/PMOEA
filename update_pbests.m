%% Update personal bests of particles in the population
function pbest=update_pbests(index,pbest,solution,numVar,namda)

  
numObjectives=2;
sum = 0;
counter = 0;
better=1;
for j=1:numObjectives
    if solution(numVar+j) <= pbest(numVar+j)
        sum=sum+1;
    end
    if solution(numVar+j) < pbest(numVar+j)
        counter=counter+1;
    end
end
if sum ==  numObjectives % /* current pop dominates pbest */<,<=
    better = 0;
else
    if sum == 0  % /* pbest dominates current pop */ >,>
        better = 1;
    elseif counter == 1  % œ‡ª•∑«÷ß≈‰ < , > sum=1
        temp1 = namda(index,1)*solution(numVar+1)+namda(index,2)*solution(numVar+2);
        temp2 = namda(index,1)*pbest(numVar+1)+namda(index,2)*pbest(numVar+2);
        if temp1<temp2
            better = 0;
        end
    end
end
if better == 0
    pbest = solution;
end

end