function Handle = P_draw(FigureData)
% ���Ƴ�ָ�����ݵ�ͼ��
% ����: FigureData,   �����Ƶ����ݾ���(��������꼯��)
%       FigureFormat, ͼ�εĸ�ʽ, ��ʡ��
% ���: Handle, ͼ�εľ��

    switch size(FigureData,2)
        case 2
            hold on;box on
            set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 13)
            if nargin < 2 || ~ischar(FigureFormat)
                Handle = plot(FigureData(:,1),FigureData(:,2),'ok','MarkerSize',6,'Marker','o','Markerfacecolor',[0 0 0]+0.7,'Markeredgecolor',[0 0 0]+0.4);
            else
                Handle = plot(FigureData(:,1),FigureData(:,2),FigureFormat);
            end
            axis tight
            xlabel('f_1')
            ylabel('f_2')
        case 3
            hold on;box on
            set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 13)
            if nargin < 2 || ~ischar(FigureFormat)
                Handle = plot3(FigureData(:,1),FigureData(:,2),FigureData(:,3),'ok','MarkerSize',8,'Marker','o','Markerfacecolor',[0 0 0]+0.7,'Markeredgecolor',[0 0 0]+0.4);
            else
                Handle = plot3(FigureData(:,1),FigureData(:,2),FigureData(:,3),FigureFormat);
            end
            axis tight
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_3')
            view(135,30)
        otherwise
            fprintf('�������2ά��3άͼ��,�޷�����.\n')
            Handle = NaN;
    end
end

