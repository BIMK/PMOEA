function Draw(population,numVar)
ParetoFront=unique(population,'rows');
M=ParetoFront(:,numVar+1);
N=ParetoFront(:,numVar+2);

Pareo_1=sort(M);
pareo_1_max=Pareo_1(size(M,1),1);
pareo_1_min=Pareo_1(1,1);
pareo_2=sort(N);
pareo_2_max=pareo_2(size(N,1),1);
pareo_2_min=pareo_2(1);
for j=1:size(M,1)
    ML(j,1)=(ParetoFront(j,numVar+1)-pareo_1_min)/(pareo_1_max-pareo_1_min);
    ML(j,2)=(ParetoFront(j,numVar+2)-pareo_2_min)/(pareo_2_max-pareo_2_min);
    if ML(j,1) == 0&&ML(j,2)==0
        ParetoFront(j,numVar+1)
        ParetoFront(j,numVar+2)
    end
end
%% save_results
P_draw(ML);