function com_max=multi_label(neighbors,a)
b=[];
for i=1:length(neighbors)
    b(i)=a(neighbors(i));
end
[~,~,e]=mode(b);
E=e{1};
com_max=E(randi(length(E)));
    


% if length(E)>1
%     com_max=a(k);
% else
%     com_max=E(randi(length(E)));
% end
% end