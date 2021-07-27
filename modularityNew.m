function Q = modularityNew(A,S,degree,edge_num)
% A = [0,1,1;1,0,0;1,0,0]; %邻接矩阵定义，与上面的例子是一致的
% S = [0,1;1,0;0,1]; %label定义
% m = sum(sum(A))/2;
% k = sum(A,2);
m = edge_num;
k = degree;
B = A - (repmat(k,[1,size(A,1)]) .* repmat(k',[size(A,1),1])) / (2*m);
Q = 1/(2*m) .* trace(S'*B*S);
end