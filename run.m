clc;
clear ;

global total_algo_time time_MOEA time_EA

k=0;           %k控制不同的网络

overlapping =0; %是否是重叠网络
% if isempty(gcp('nocreate'))
%     parpool(28);
% end
switch k

%6个真实网络;
case 0
  isRealWorld = 1;
Date={'karate','dolphin','football','netscience_remove','CQC','blogs','hepth','hepth1','CA-AstroPh_18772_396160','CA-CondMat_23133_186936','Brightkite_58228','196591'};
    for i=12:12
        path = sprintf('RealWorld/%s.txt',Date{i});
        name=Date{i};
        real_path=sprintf('RealWorld/real_label_%s.txt',Date{i});
          for c=1:1
          total_algo_time = 0;
%           tic;
              
%          PMOEA_MOEA(path,name,real_path,overlapping,isRealWorld,1,-1,c);
          PMOEA_EA(path,name,real_path,overlapping,isRealWorld,c);
          
%           time_used_total = toc
          filename_time = sprintf('result_%s_%d.txt', name, c);
		  savedata1(filename_time, [total_algo_time time_MOEA time_EA]);
          end
     end
    
    
 
 %% 实验1000_20_50_20_100_0.1_0.8
     case 1000
     Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
     isRealWorld = 0;
     overlapping = 0;
     for i=1:1
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/1000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/1000/%s/community.dat',Date{i}); 
       name=sprintf('1000_%s',Date{i});
       for c=1:1
           PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
           PMOEA_EA(path,name,real_path,0,0,c);
       end
     end
     
     %% 实验1000_20_50_20_100_0.1_0.8
     case 10001
     Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
     isRealWorld = 0;
     overlapping = 0;
     for i=2:8
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/1000_test/%s/network.dat',Date{i});
       real_path=sprintf('LFR/1000_test/%s/community.dat',Date{i}); 
       name=sprintf('1000_test_%s',Date{i});
       for c=1:2
%            PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
           PMOEA_EA(path,name,real_path,0,0,c);
       end
     end
    
 %% 实验2000_20_50_20_100_0.1_0.8
     case 2000
         
     Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
     isRealWorld = 0;
     overlapping = 0;
    for i=7 :8
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/2000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/2000/%s/community.dat',Date{i}); 
       name=sprintf('2000_%s',Date{i});
       for c=1:1
            PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
            PMOEA_EA(path,name,real_path,0,0,c);
       end
    end
    
    %% 实验5000_20_50_20_100_0.1_0.8
     case 5000
         
     Date={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'};
     isRealWorld = 0;
     overlapping = 0;
    for i=6:6
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/5000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/5000/%s/community.dat',Date{i}); 
       name=sprintf('5000_%s',Date{i});
       t1 = clock;

       for c=1:1
             PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
            PMOEA_EA(path,name,real_path,0,0,c);
       end
       t2 = clock;
 etime(t2,t1)

    end
    
    
    case 10000

    Date={'0.1','0.25','0.3','0.4','0.5','0.6','0.7','0.8'};
    isRealWorld = 0;
    overlapping = 0;
    for i=2 :2
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/10000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/10000/%s/community.dat',Date{i});
       name=sprintf('10000_%s',Date{i});
       for c=1:1
%             PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
            PMOEA_EA(path,name,real_path,0,0,c);
       end
    end
    
    case 20000

    Date={'0.1','0.25','0.3','0.4','0.5','0.6','0.7','0.8'};
    isRealWorld = 0;
    overlapping = 0;
    for i=2 :2
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/20000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/20000/%s/community.dat',Date{i});
       name=sprintf('20000_%s',Date{i});
%       PMOEA_MOEA(path,name,real_path,0,0,201,400);
      PMOEA_EA(path,name,real_path,0,0);
    end

    case 30000

    Date={'0.1','0.25','0.3','0.4','0.5','0.6','0.7','0.8'};
    isRealWorld = 0;
    overlapping = 0;
    for i=2 :2
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/30000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/30000/%s/community.dat',Date{i});
       name=sprintf('30000_%s',Date{i});
      PMOEA_MOEA(path,name,real_path,0,0,116,250);
%       PMOEA_EA(path,name,real_path,0,0);
    end
    
    case 40000

    Date={'0.1','0.25','0.3','0.4','0.5','0.6','0.7','0.8'};
    isRealWorld = 0;
    overlapping = 0;
    for i=2 :2
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/40000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/40000/%s/community.dat',Date{i});
       name=sprintf('40000_%s',Date{i});
%       PMOEA_MOEA(path,name,real_path,0,0,401,500);
      PMOEA_EA(path,name,real_path,0,0);
    end
    
    
    case 50000

    Date={'0.1','0.25','0.3','0.4','0.5','0.6','0.7','0.8'};
    isRealWorld = 0;
    overlapping = 0;
    for i=2 :2
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR/50000/%s/network.dat',Date{i});
       real_path=sprintf('LFR/50000/%s/community.dat',Date{i});
       name=sprintf('50000_%s',Date{i});
       for c=1:5
            PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
            PMOEA_EA(path,name,real_path,0,0,c);
       end
    end
    
            case 10005000

    Date={'1000','2000','3000','4000','5000'};
    isRealWorld = 0;
    overlapping = 0;
    for i=1:5
       clearvars -EXCEPT Date i change Noture edge_add;
       path=sprintf('LFR_1000-5000/%s.txt',Date{i});
       real_path=sprintf('LFR_1000-5000/real_%s.txt',Date{i});
       name=sprintf('LFR_1000-5000_%s',Date{i});
       for c=1:5
            PMOEA_MOEA(path,name,real_path,0,0,1,-1,c);
            PMOEA_EA(path,name,real_path,0,0,c);
       end
    end
    
end

quit
 
