function f = evaluate_objective(x,AdjMatrix,ll)
%%  ll=sum(sum(AdjMatrix,1))
M=1;                       %%M=1��ʾΪ��ǩ������M=0��ʾΪ�ڽӱ���
f=[];
%x= [4,2,4,4,4,4,4,4,2,10,4,4,4,4,15,34,17,4,34,4,34,4,34,33,32,33,27,33,29,27,2,4,2,2];
if  M==1
    clu_assignment=x(1:size(AdjMatrix,1));%��ǩ��������
end
if M==0
    clu_assignment = decode(x(1:size(AdjMatrix,1)));%�ڽӱ������
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���߱�ʾ�Ļ���ת��ɵ�Ļ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f(1) = community_fitness_mod2(AdjMatrix,clu_assignment,ll);
f(2) = community_score(AdjMatrix,clu_assignment,ll);

end