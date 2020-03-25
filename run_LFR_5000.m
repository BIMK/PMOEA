function  [overlapping]=run_LFR_5000(fid)

     Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
     isRealWorld = 0;
     overlapping = 0;
     
    for i=1:length(Date)
        
       clearvars -EXCEPT Date i change Noture edge_add isRealWorld overlapping fid;
       path=sprintf('LFR/5000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/5000/%s/community.dat',Date{i}); 
       name=sprintf('5000_%s',Date{i});
       
       save_root = ['result_statistics/',name];
      if ~isdir(save_root) %判断路径是否存在
           mkdir(save_root);
      end
      
      t_start = clock;
fprintf(fid, '当前时间是%s\n',datestr(now));       
     
       for c=1:1
           
           
           t0 = clock;
           
         
           PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
           t1 = clock;
           
           PMOEA_EA(path,name,real_path,0,0,c);
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
            
            STR = ['Network:',name,', run:',num2str(c),'/20',', this run used time:',num2str(etime(clock,t0)),', total used time:',num2str(etime(clock,t_start))];
            fprintf(fid,'%s\n',STR);
            
            
       end
     

    end
    

fclose(fid);
end
