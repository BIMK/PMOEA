function chromosomes=update_neighbour(idealp,chromosomes,child,i_neighbour_index,weight,niche) 
% niche=length(i_neighbour_index);
for i=1:niche
    neighbours = i_neighbour_index(i);
    f1 = scalar_func(chromosomes{i}(1,end-1:end),idealp,weight(neighbours,:));
    f2 = scalar_func(child{1}(1,end-1:end),idealp,weight(neighbours,:));
  if f2 < f1    %¸üÐÂÁÚÓò½â
      chromosomes(neighbours) = child;
  end
end