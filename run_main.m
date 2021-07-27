function [P] = run_main(N)

    fid=fopen('journal','a');
fprintf(fid, '############开始运行############%s\n',datestr(now));
    [a] = karate(fid);
global modularity
modularity = @modularity_;
   %[overlapping]=run_LFR_5000(fid);
    
    for N=1:4
    %    [N_nodes]=run_LFR_5W_2_20W(N,fid);
    end    
    
    fclose(fid);
    
    P=0;
    


end
