%==========================================================
% MATLAB job submission script: batch2.m
%==========================================================
%% Set Workers
cores = 4;    %设置作业请求的核数
workers = cores-1;    %设置matlab提交作业时的workers数

%% Import Cluster Profile
JobPath = pwd    %获取当前作业路径
allName = parallel.clusterProfiles();    %获取当前matlab的集群配置文件清单
if ~ismember('Slurm_ahuhpc',allName)    %判断是否存在名称为Slurm_ahuhpc的配置文件
parallel.importProfile([JobPath,'/Slurm_ahuhpc'])    %如果不存在则导入作业路径下的Slurm_ahuhpc配置文件
end

%% Culster Settings
tic;
c = parcluster('Slurm_ahuhpc');    %调用Slurm_ahuhpc配置文件创建集群并行池
c.ClusterMatlabRoot = '/Share/apps/matlab/R2019b';    %设置matlab的安装路径
c.SubmitArguments = '--ntasks-per-node=28'   %添加调度系统作业提交参数,--ntasks-per-node参数限制每节点的任务数
fprintf('Creating the cluster: %.2f sec\n', toc);

%% Submit Job
tic;
j = batch(c, @sum2, 1, {10}, 'pool', workers)    %batch参数，1、并行池 2、提计算的函数 3、函数的输出参数个数 4、函数的输入参数 5、并行池参数 6、子任务数
fprintf('Submitting the job: %.2f sec\n', toc);

tic;
j.wait   %等待所有并行计算的子任务完成
fprintf('Waiting time: %.2f sec\n', toc);

fprintf('\nJob result output\n');
j.fetchOutputs{:}    %获取作业中函数的输出结果(需要函数定义输出参数)

%% Print Job info
fprintf('\nJob Properties\n');
fprintf('Create Time: %s \n', j.CreateTime);
fprintf('Submit Time: %s \n', j.SubmitTime);
fprintf('Start Time: %s \n', j.StartTime);
fprintf('Finish Time: %s \n', j.FinishTime);
fprintf('\nEnd\n');

exit

