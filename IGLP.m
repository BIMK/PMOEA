function gene = IGLP(numVar,AdjMatrix,U)
degree=single(sum(AdjMatrix,1));
if U==0     %%±êÇ©´«²¥
num=1;
gene = single(1:numVar);
a=randperm(numVar);
for n = 1 : 3    
    for  i = 1 : numVar    
        %%%%%%NeighborSize = degree(i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         NeighborSize=degree(a(i));
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if NeighborSize == 0   
            %			gene(i) = 0;            
        else            
            if NeighborSize == 1   
                neighbours=find(AdjMatrix(a(i),:)==1);
                gene(a(i)) = gene(neighbours(1));                
            else  
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 neighbours=find(AdjMatrix(a(i),:)==1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 sum = 0;
                maxr = -1;%//record index of i's neighbour which ...
                label = -1;
                temp = 1;                
                for j = 1 : NeighborSize                   
                    counter = 1; %//record no. of nodes that has same label with j 
                    for k = j + 1 : NeighborSize 
                        p = gene(neighbours(j));
                        q = gene(neighbours(k));
                        
                        if p == q
                            counter = counter + 1;
                        end
                    end %//end k                    
                    if temp < counter                       
                        maxr = j;
                        temp = counter;
                    end
                end %//end j                
                for l =1 : NeighborSize                    
                    u = gene(neighbours(1));
                    v = gene(i);
                    if  u == v                        
                        label = u;
                    end
                end %//end l
                if label ~= -1 && maxr == -1                    
                    gene(a(i)) = label;                   
                else                    
                    if maxr ~= -1                        
                        gene(a(i)) = gene(neighbours(maxr)); 
                        if(length(neighbours)>=3) 
                            gene(a(i))=a(i); 
                        end
                    else                       
                        randneighbor = randi(NeighborSize);
                        %randneighbor = 2;%test
                        gene(a(i)) = gene(neighbours(randneighbor));  
                        if(length(neighbours)>=3) 
                            gene(a(i))=a(i); 
                        end
                    end
                end
            end % if NeighborSize == 1 
        end % if NeighborSize == 0
    end %//end i
end %//end n
else 
    for i=1:numVar
     gene(i)=single(i);
    end
end