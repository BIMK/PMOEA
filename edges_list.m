function [et] = edges_list(A,V)
%EDGELIST Generate the edges list for each node.
et = [];
for i = 1:V
    et(i).e = [];
    et(i).n = 0;
    for j = 1:V
        if (A(i,j) == 1)
            et(i).e = [et(i).e j];
            et(i).n = et(i).n + 1;
        end
    end   
end
end