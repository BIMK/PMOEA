function [a] = karate(fid)
    isRealWorld = 1;
    global total_algo_time time_MOEA time_EA
    overlapping = 0;
 Date = {'karate', 'dolphin', 'football', 'netscience_remove', 'CQC', 'blogs', 'hepth', 'hepth1', 'CA-AstroPh_18772_396160', 'CA-CondMat_23133_186936', 'Brightkite_58228', '196591'};
    %Date={'karate','dolphin','football','netscience_remove','CQC','blogs','Hepth','hepth1','CA-AstroPh','CA-CondMat','Brightkite_edges','Gowalla_edges'};
    for i=1:11
        
      global modularity
      if i==12 %最大的数据集需要特殊处理
  	modularity = @modularity_revised;
      else
 	modularity = @modularity_;
      end
      save_root = ['result_statistics/',Date{i}];
      if ~isdir(save_root) %判断路径是否存在
           mkdir(save_root);
      end     
      
      path = sprintf('RealWorld/%s.txt',Date{i});
      name=Date{i};
      real_path=sprintf('RealWorldreal_label_%s.txt',Date{i});
      
       t_start = clock;
       
      
        for c=1:1
            total_algo_time = 0;
       
            
            t0 = clock;            
            
            PMOEA_MOEA(path,name,real_path,overlapping,isRealWorld,1,-1,c);
            %break
            
            t1 = clock;
            PMOEA_EA(path,name,real_path,overlapping,isRealWorld,c);
            t_end = clock;
            
            
            third_used_time = etime(t_end,t1);
            total_used_time = etime(t_end,t0);
            
            save_root = ['result_statistics/',name,'/third_used_time'];
            if ~isdir(save_root) %判断路径是否存在
                mkdir(save_root);
            end
            dlmwrite([save_root,'/third_used_time_',num2str(c),'.txt'],third_used_time,' ');
            
            save_root = ['result_statistics/',name,'/total_used_time'];
            if ~isdir(save_root) %判断路径是否存在
                mkdir(save_root);
            end
            dlmwrite([save_root,'/total_used_time',num2str(c),'_.txt'],total_used_time,' ');
            
            %           time_used_total = toc
            filename_time = sprintf('result_%s_%d.txt', name, c);
            %savedata1(filename_time, [total_algo_time time_MOEA time_EA]);
            a= [total_algo_time time_MOEA time_EA];
            
            STR = ['Network:',name,', run:',num2str(c),'/20',', this run used time:',num2str(etime(clock,t0)),', total used time:',num2str(etime(clock,t_start))];
            fprintf(fid,'%s\n',STR);
            
            
        end
   end
  

end
