function chromosomes=update_problem(idealp,chromosomes,child,i_neighbour_index,index,weight,numVar) 


niche=length(i_neighbour_index);
for i=1:niche
    neighbours=i_neighbour_index(i);
    f1=scalar_func(chromosomes(neighbours,numVar+1:numVar+2),idealp,weight(neighbours,:));
    f2=scalar_func(child(numVar+1:numVar+2),idealp,weight(neighbours,:));
  if f2<f1    %¸üÐÂÁÚÓò½â
      chromosomes(neighbours,:)=child;
  end
end
      