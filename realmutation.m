%% NBM mutation used in the paper
function f = realmutation(ind,mutate_posibility,AdjMatrix,similarity)
%%ÿ�����Ӵ���һ�����绮�֣�ÿ�����Ӷ���N�����㣨����λ������Ӧ����ڵ�������ɣ�����λֵ��ʾ���绮�֣�
%%�������Ӷ�Ӧͬһ���ڽӾ����������˽ṹ��
numVar=size(AdjMatrix,1);
degree=sum(AdjMatrix,1);
for j =1:numVar     %% ��ind������ӵ�ÿ������λ�������j���ڵ㣩������
                    %%ÿ����������0��1֮���α���������������С��ͻ�����pm��
                    %%��NBM����Ӧ���ڸû��򣬼��� �����ǩ��ʶ��������������ھӡ�
% %     neighbours=find(AdjMatrix(j,:)==1);%��j��������ھӼ���
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�и���ֻѡȡ���ƶȱȽϴ���ھ�
%     if rand <= mutate_posibility
        neighbours=find(similarity(j,:)>=0.10);
%     else
%         neighbours=find(similarity(j,:)>=0);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rand <= mutate_posibility
%     if rnd_uni() <= mutate_posibility
        indentifier = ind(j);
        %         negative_n = node(j).neighbours_n.size();
% %         for i = 1:degree(j)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for i = 1:size(neighbours,2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            neighborX = neighbours(i);
            %%
            %%���ھӵ����ƶ�С����ֵ���򲻴���
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if similarity(j,neighborX)>=0.3
                ind(neighborX) = indentifier;
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %ind(neighborX) = indentifier;
        end
    end
end
f=ind;
end