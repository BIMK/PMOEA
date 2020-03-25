function Q = modularity_revised(solution,AdjacentMatrix,degree,edge_num )
%% 获取无向无权网络的某种划分的模块

m=max(solution);
edge_num = edge_num*2;
Q=0;
for i=1:m
    index = single(find(solution==i));
    if ~isempty(index)

        d = degree(index);
        e = d';              
	if length(d)<1000
           temp_sum=  sum(sum(AdjacentMatrix(index,index)-d*e/edge_num,1));
        else
        temp_sum=Q_sum(AdjacentMatrix(index,index),d,edge_num);
        end
        Q=Q+temp_sum;
    end
end

clear index solution d e;
Q=Q/edge_num;

end


%function temp_sum=Q_sum(AdjacentMatrix_index,d,edge_num)
%
%temp_sum = 0;
%for i=1:length(d)
%    D_matrix_i_row = (d(i).*d')./edge_num;    
   % D_matrix_i_row = (d(i)*d')/edge_num;    
%    temp_sum = temp_sum+ sum(AdjacentMatrix_index(i,:)-D_matrix_i_row);   
%end

%end
function temp_sum=Q_sum(AdjacentMatrix_index,d,edge_num)

temp_sum = 0;

step=round(length(d)/10);
ternal = [0:step:length(d)];
if ternal(end)~=length(d)
    ternal = [ternal length(d)];    
end
total_step = length(ternal);

for j=1:total_step-1 
    i = [ternal(j)+1:ternal(j+1)];
    D_matrix_i_row = (d(i)*d')./edge_num;    
    temp_sum = temp_sum+ sum(sum(AdjacentMatrix_index(i,:)-D_matrix_i_row));   
end
end
