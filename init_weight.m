function [weights,neighbours] = init_weight(popsize, niche)
% init_weights and neighbors.

weights = [];
for i=0:popsize-1
    weight=zeros(1,2);
    weight(1)=i/(popsize-1);
    weight(2)=(popsize-i-1)/(popsize-1);
    weights = [weights;weight];
end
weights=single(weights);

%Set up the neighbourhood.
leng=size(weights,1);
distanceMatrix=zeros(leng, leng);
for i=1:leng
    for j=i+1:leng
        A=weights(i,:)';
        B=weights(j,:)';
        distanceMatrix(i,j)=(A-B)'*(A-B);
        distanceMatrix(j,i)=distanceMatrix(i,j);
    end
    [s,sindex]=sort(distanceMatrix(i,:));
    neighbours(i,:)=sindex(1:niche);
end
   neighbours=single(neighbours);
end