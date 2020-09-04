"""Methods for scheduling XGAP through SGE"""

from subprocess import check_output, CalledProcessError

def submit_job(job_name, job, mem, runtime, log_output, num_tasks=1, highp=False,
               queue=None, hold_id=None):
  """Submits job with given parameters to SGE queue"""
  cmd = ["qsub",
         "-N", job_name,
         "-cwd", "-V",
         "-j", "y",
         "-o", log_output]
  resource_str = "h_data={},h_rt={}".format(mem, runtime)
  if highp:
    resource_str += ",highp"
  cmd.append("-l")
  cmd.append(resource_str)
  if queue:
    cmd.append("-q")
    cmd.append(queue)
  if num_tasks > 1:
    cmd.append("-t")
    cmd.append("1-{}:1".format(num_tasks))
  if hold_id:
    cmd.append("-hold_jid")
    cmd.append(hold_id)
  cmd = cmd + job.split()
  try:
    output = check_output(cmd)
  except CalledProcessError:
    raise Exception("Could not submit job")
  try:
    output_str = output.decode("utf-8")
  except AttributeError:
    output_str = output
  print(output_str)
  # Assumes last line of stdout is 
  # "Your [job] [job_id] ([job_name]) has been submitted"
  job_id = output_str.strip().split()[-5].split('.')[0]
  if not job_id.isnumeric():
    raise Exception("Job ID not found")
  return job_id

def delete_job(job_id):
  cmd = ["qdel", job_id]
  try:
    output_str = check_output(cmd)
  except CalledProcessError:
    raise Exception("Could not delete job {}".format(job_id))
  print(output_str)

if __name__ == "__main__":
  pass
  #delete_job(submit_job("test_job", "test.sh", "10G", "12:00:00", "log.txt"))

