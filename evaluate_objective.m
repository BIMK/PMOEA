function f = evaluate_objective(x,AdjMatrix,ll)
%%  ll=sum(sum(AdjMatrix,1))
M=1;                       %%M=1表示为标签传播，M=0表示为邻接编码
f=[];
%x= [4,2,4,4,4,4,4,4,2,10,4,4,4,4,15,34,17,4,34,4,34,4,34,33,32,33,27,33,29,27,2,4,2,2];
if  M==1
    clu_assignment=x(1:size(AdjMatrix,1));%标签传播解码
end
if M==0
    clu_assignment = decode(x(1:size(AdjMatrix,1)));%邻接编码解码
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%将边表示的划分转变成点的划分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f(1) = community_fitness_mod2(AdjMatrix,clu_assignment,ll);
f(2) = community_score(AdjMatrix,clu_assignment,ll);

end