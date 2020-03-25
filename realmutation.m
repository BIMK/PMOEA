%% NBM mutation used in the paper
function f = realmutation(ind,mutate_posibility,AdjMatrix,similarity)
%%每个粒子代表一个网络划分，每个粒子都由N个顶点（基因位）（对应网络节点数）组成，基因位值表示网络划分，
%%所有粒子对应同一个邻接矩阵（网络拓扑结构）
numVar=size(AdjMatrix,1);
degree=sum(AdjMatrix,1);
for j =1:numVar     %% 对ind这个粒子的每个基因位（网络第j个节点）操作，
                    %%每个基因生成0和1之间的伪随机数，如果随机数小于突变概率pm，
                    %%则将NBM过程应用于该基因，即， 将其标签标识符分配给其所有邻居。
% %     neighbours=find(AdjMatrix(j,:)==1);%第j个顶点的邻居集合
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %有概率只选取相似度比较大的邻居
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
            %%若邻居的相似度小于阈值，则不传播
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