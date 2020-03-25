function [ Labelnew ] = LPA( adjacent_matrix,label )
    if nargin<2
        label = 1:size(adjacent_matrix,2);
    end
    N = size(adjacent_matrix,2);
   
    Label1 = label;
    Label2 = Label1;
    Labelnew = Label1;
    flag=1;
    while(1)
        for i=1:N
            nb_lables = Labelnew(find(adjacent_matrix(i,:)==1));%�ҵ��ھ��±��Ӧ�ı�ǩ
            if size(nb_lables,2)>0
                x = tabulate(nb_lables);
                max_nb_labels = x(find(x(:,2)==max(x(:,2))),1);
                Labelnew(i) = max_nb_labels(randi(length(max_nb_labels))); 
            end
        end
        % ��������,Ԥ����Ծ
        if all(Labelnew==Label1)||all(Labelnew==Label2)
            break;
        else
            if flag==1
                Label1 = Labelnew;
                flag=0;
            else
                Label2 = Labelnew;
                flag=1;
            end
        end
    end
    
end

