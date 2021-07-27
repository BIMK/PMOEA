function [] = subjob()
  workers = 4-1;

  tic;
  c = parcluster('Slurm_ahuhpc');
  c.ClusterMatlabRoot = '/Share/apps/matlab/R2019b';
  c.SubmitArguments = ['--ntasks-per-node=2 -e slurm-%j.out -o slurm-%j.out -J ',getenv('SLURM_JOB_NAME'),'-sub']
  fprintf('Creating the cluster: %.2f sec\n', toc);

  tic;
  j = batch(c, 'sum2', 1, {100}, 'pool', workers)
  fprintf('Submitting the job: %.2f sec\n', toc);

  tic;
  j.wait
  fprintf('Waiting time: %.2f sec\n', toc);

  fprintf('\nJob result output\n');
  j.fetchOutputs{:}

  fprintf('\nJob Properties\n');
  fprintf('Create Time: %s \n', j.CreateTime);
  fprintf('Submit Time: %s \n', j.SubmitTime);
  fprintf('Start Time: %s \n', j.StartTime);
  fprintf('Finish Time: %s \n', j.FinishTime);
  fprintf('\nEnd\n');

end
