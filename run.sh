#!/bin/bash
if [ $# -lt 5 ]; then
	echo -e "错误:  参数不足"
	echo -e "请按照 [脚本] [m文件名] [返回值个数] [输入参数] [核数] [单节点核数限制] <作业名> 这种格式添加参数来运行。"
    exit 1
else
	echo "输入的m文件:    " $1
	echo "返回值个数:     " $2
	echo "输入参数:       " $3
	echo "申请的核数:     " $4
	echo "单节点核数限制: " $5
fi

if [ ! "$6" == "" ]; then
	JOBNAME="#SBATCH -J $6"
	echo "作业名:         " $6
	echo "子作业名:       " $6"-sub"
fi

cat > run.sbatch <<EOF
#!/bin/bash
#SBATCH -e slurm-%j.out
#SBATCH -o slurm-%j.out
#SBATCH -p COMPUTE
$JOBNAME

module load /Share/apps/matlab/R2019b
/Share/apps/matlab/R2019b/bin/matlab -nodesktop -nosplash -nodisplay -r "subjob;exit"
EOF

cat > subjob.m <<EOF
function [] = subjob()
  workers = $4-1;

  tic;
  c = parcluster('Slurm_ahuhpc');
  c.ClusterMatlabRoot = '/Share/apps/matlab/R2019b';
  c.SubmitArguments = ['--ntasks-per-node=$5 -e slurm-%j.out -o slurm-%j.out -J ',getenv('SLURM_JOB_NAME'),'-sub']
  fprintf('Creating the cluster: %.2f sec\\n', toc);

  tic;
  j = batch(c, '$1', $2, {$3}, 'pool', workers)
  fprintf('Submitting the job: %.2f sec\\n', toc);

  tic;
  j.wait
  fprintf('Waiting time: %.2f sec\\n', toc);

  fprintf('\\nJob result output\\n');
  j.fetchOutputs{:}

  fprintf('\\nJob Properties\\n');
  fprintf('Create Time: %s \\n', j.CreateTime);
  fprintf('Submit Time: %s \\n', j.SubmitTime);
  fprintf('Start Time: %s \\n', j.StartTime);
  fprintf('Finish Time: %s \\n', j.FinishTime);
  fprintf('\\nEnd\\n');

end
EOF
 
sbatch run.sbatch

