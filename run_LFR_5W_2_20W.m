function [N_nodes] = run_LFR_5W_2_20W(N, fid)

    number_nodes = [1 2 3 4] .* 50000;
    N_nodes = number_nodes(N);

    Date = {'0.1', '0.25', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8'};
    isRealWorld = 0;
    overlapping = 0;

    for i = 2:2

        clearvars -EXCEPT Date i change Noture edge_add N_nodes fid;
        path = sprintf('LFR/%d/%s/network.dat', N_nodes, Date{i});
        real_path = sprintf('LFR/%d/%s/community.dat', N_nodes, Date{i});
        name = sprintf('%d_%s', N_nodes, Date{i});

        save_root = ['result_statistics/', name];

        if ~isdir(save_root)%判断路径是否存在
            mkdir(save_root);
        end

        t_start = clock;
        fprintf(fid, '当前时间是%s\n', datestr(now));

        for c = 1:20

            t0 = clock;
            PMOEA_MOEA(path, name, real_path, 0, 0, 1, -1, c);
            t1 = clock;

            PMOEA_EA(path, name, real_path, 0, 0, c);
            t_end = clock;

            third_used_time = etime(t_end, t1);
            total_used_time = etime(t_end, t0);

            save_root = ['result_statistics/', name, '/third_used_time'];

            if ~isdir(save_root)%判断路径是否存在
                mkdir(save_root);
            end

            dlmwrite([save_root, '/third_used_time_', num2str(c), '.txt'], third_used_time, ' ');

            save_root = ['result_statistics/', name, '/total_used_time'];

            if ~isdir(save_root)%判断路径是否存在
                mkdir(save_root);
            end

            dlmwrite([save_root, '/total_used_time', num2str(c), '_.txt'], total_used_time, ' ');

            STR = ['Network:', name, ', run:', num2str(c), '/20', ', this run used time:', num2str(etime(clock, t0)), ', total used time:', num2str(etime(clock, t_start))];
            fprintf(fid, '%s\n', STR);

        end

    end

    fclose(fid);
end
