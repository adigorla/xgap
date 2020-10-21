"""Methods for scheduling XGAP through SLURM"""

from subprocess import check_output, CalledProcessError
from sys import stdout

### currently only compatible withe SLURM 18.08.6-2 and pcluster 2.4.1

def submit_job(job_name, job, mem, runtime, log_output, num_tasks=1, highp=False,
               queue=None, hold_ids=None):
  """Submits job with given parameters to SLURM queue"""
  if log_output != "/dev/null":
    log_output = "{}/%x.%j.log".format(log_output)
  cmd = ["sbatch",
         "--job-name={}".format(job_name),
         "--export=ALL",
         "--output={}".format(log_output),
         "-t", runtime,
        "-N", "1",
        "-n", "1",
#        "-c", 1,
        "--mem-per-cpu={}".format(mem)]
  if highp:
    pass
  if queue:
    cmd.append("-p")
    cmd.append(queue)
  if num_tasks > 1:
    cmd.append("--array=1-{}:1".format(num_tasks))
  if hold_id:
    hold_ids_str = ":".join(ids.strip() for ids in hold_ids)
    cmd.append("-d")
    cmd.append("afterany:{}".format(hold_ids_str))
  cmd = cmd + ["--wrap=/usr/bin/python3 {}".format(job)]
  try:
    output = check_output(cmd)
  except CalledProcessError:
    raise Exception("Could not submit job")
  # Assumes stdout upon sucessful submission is 
  # "Submitted batch job [SLURM_JOBID]"
  try:
    output_str = output.decode("utf-8")
  except AttributeError:
    output_str = output
  stdout.write(output_str)
  stdout.flush()
  #CHANGE THIS LINE( V ) IF "Job ID not found" exception thrown
  job_id = output_str.strip().split()[-1]
  if not job_id.isnumeric():
    raise Exception("Job ID not found, conisder editing line 37 in ~/xgap/scheduler/slurm.py")
  return job_id

def delete_job(job_id):
  cmd = ["scancel", job_id]
  try:
    output_str = check_output(cmd)
  except CalledProcessError:
    raise Exception("Could not delete job {}".format(job_id))
  print(output_str)

def dummy():
    print("in function name space!")

if __name__ == "__main__":
  pass
  #delete_job(submit_job("test_job", "test.sh", "10G", "12:00:00", "log.txt"))

