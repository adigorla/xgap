"""Methods for checkpoint system in pipeline"""

from abc import ABCMeta, abstractmethod, abstractproperty
from os import path
from os.path import exists, getsize
import os
from subprocess import Popen, PIPE, check_output, CalledProcessError
from sys import argv, stdout
import smtplib, ssl
from datetime import datetime

__all__ = [ "Task" ]

# How many times a task can start before being aborted
start_limit = 2

# Max number of failed tasks that will be rerun.
# If more, analysis will be aborted
rerun_limit = 10

class Task(object):
  __metaclass__ = ABCMeta

  def __init__(self, sample_id, sample_path, config, config_path):
    self.sample_id = sample_id
    self.sample_path = sample_path
    self.config = config
    self.config_path = config_path
    self.xgap_path = self.config["xgap_path"]
    self.log_dir = "{}/xgap_logs/{}".format(config["log-dir"],self.sample_id)
    self.out_dir = "{}/xgap_output/{}".format(config["output-dir"],self.sample_id)
    self.output_paths = self.gen_output_paths()
    self.log_paths = self.gen_log_paths()
  
  #TODO: make dedicated error handlers for each module.
  def check(self):
    """Return true if task is complete"""
    return check_step(self.output_paths, self.log_paths,
                      self.success_terms, self.error_terms)
  @abstractproperty
  def error_terms(self):
    """Define abstract member error_terms"""
    pass

  @abstractproperty
  def success_terms(self):
    """Define abstract member success_terms"""
    pass

  @abstractmethod
  def gen_output_paths(self):
    """Initialize abstract member output_paths (list of strings)"""
    pass

  @abstractmethod
  def gen_log_paths(self):
    """Initialize abstract member log_paths (list of strings)"""
    pass

  @abstractmethod
  def run(self):
    """Perform task"""
    pass

  @abstractmethod
  def clean_up(self):
    """Method run after run() successfully completes"""
    pass
  @abstractmethod
  def mail_notify(self, msg):
    """Method to send email noification of msg
       Executed when program is aborted or pipeline run is complete,
       Also deletes file with email config info
       NOTE: only compatable with gmail SEND FROM addresses for now
    """
    mailf = self.log_dir+"/mail.dat"
    sbj = str(self.sample_id) + " Analysis " + msg[0]
    bdy = []
    time = datetime.now()
    bdy.append(sbj+" at "+time.strftime("%Y-%m-%d %H:%M:%S")+" local time\n")
    if (msg[0]is not None) and (msg[0].lower() == "aborted"):
      bdy.append(sbj+" during "+msg[1]+" \n")
      bdy.append(msg[2])
      bdy.append("\nThe following region indices failed:")
      if not isinstance(msg[3], list):
        msg[3] = [msg[3]]
      bdy += msg[3]
    fmsg = "Subject: {}\n\n".format(sbj)
    for itm in bdy:
      fmsg += "{}\n".format(itm)
      
    if _check_file(mailf):
      with open(mailf, 'r') as mf:
        info = [line.strip() for line in mf]

      try:
        port = 465
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
          server.login(info[1], info[3])
          server.sendmail(info[1], info[2], fmsg)
        status = "Email notification sent to :{}\n".format(info[2])

      except smtplib.SMTPAuthenticationError:
        status = "Email notification:  Failed, wrong sender email password OR mail config info is not correct"

      finally:
        os.remove(mailf)
        return status
    
    else:
      return "Email notification: Failed, mail config info not found"

def check_step(output_paths, log_paths, success_terms, error_terms):
  """Returns true if step is successful and outputs all issues if they exist

  Success entails:
    1: all files in output_paths exist and are nonempty
    2: each log in log_paths has all success_terms and no error_terms
    3: output_paths and log_paths are both nonempty

  Args:
    output_paths: List of output path lists for each task index
    log_paths: List of log path lists for each task index
    success_terms: List of strings that should be in log file
    error_terms: List of strings that indicate an error occurred
  """
  failed_indices = set()
  if len(output_paths) != len(log_paths):
    raise Exception("Output or log files not specified for all tasks")
  n_tasks = len(output_paths)
  for index in range(n_tasks):
    for task_output_path in output_paths[index]:
      if not _check_file(task_output_path):
        print("Output file at {} does not exist or is empty".format(task_output_path))
        failed_indices.add(index)
    for task_log_path in log_paths[index]:
      if not _check_file(task_log_path):
        print("Log file at {} does not exist or is empty".format(task_log_path))
        failed_indices.add(index)
        continue
      if _check_log(task_log_path, success_terms, error_terms):
        failed_indices.add(index)
  return sorted(failed_indices)
  # Check that output is valid
  """
  success = True
  if not output_paths:
    print("No output files listed")
    success = False
  if not log_paths:
    print("No log files listed")
    success = False
  for output_path in output_paths:
    if not _check_file(output_path):
      success = False
      print("Output file at {} does not exist or is empty".format(output_path))
  # Check that logs are normal
  for log_path in log_paths:
    if _check_log(log_path, success_terms, error_terms):
      success = False

  return success
  """

def _check_file(file_path):
  """Returns true if file at file_path exists and is nonempty"""
  return (os.path.exists(file_path) and os.path.getsize(file_path) > 0)
  #NEED TO FIX file check

def _check_log(log_path, success_terms, error_terms):
  """Checks log for specified success and error terms, ignoring case

  Note that terms should not span multiple lines or occur multiple times
  on a line :(

  Args:
    log_path: path to log_file
    success_terms: List of strings that should be in log file
    error_terms: List of strings that indicate an error occurred
  Returns:
    0 if log file indicates sucess
    -1 if error term found
    n if n success terms were not found
  """
  remaining_success_terms = list(success_terms)
  terms_to_remove = set()
  with open(log_path, 'r') as log_file:
    for line_num, line in enumerate(log_file):
      line = line.lower()
      for term in error_terms:
        if line.find(term.lower()) != -1:
          print("Found error in line {} of ({})".format(line_num, log_path))
          print("\t\"{}\"".format(line))
          return -1
      for term in success_terms:
        if line.find(term.lower()) != -1:
          if not term in terms_to_remove:
            terms_to_remove.add(term)
    for term in terms_to_remove:
      remaining_success_terms.remove(term)
  if remaining_success_terms:
    print("\nFollowing terms not found in log ({}):".format(log_path))
    for term in remaining_success_terms:
      print("\t{}".format(term))
    return len(remaining_success_terms)
  return 0

def check_disk_usage(config):
  """Returns true if disk usage under threshold, false otherwise"""
  disk_id = config["output-disk-identifier"]
  disk_threshold = config["free-disk-threshold"]
  disk_usage = _disk_usage(disk_id)
  if disk_usage > disk_threshold:
    print("Disk usage ({}%) above threshold of {}%".format(disk_usage,
                                                           disk_threshold))
    return False
  return True

def _disk_usage(disk_id):
  """Returns disk usage percentage (ex. "99%" returns 99) for given disk_id

  Calls df and greps for disk_id. If other than 1 entry, raises exception.
  """
  df_ps = Popen(('df'), stdout=PIPE)
  cmd = ["grep", disk_id]
  try:
    output = check_output(cmd, stdin=df_ps.stdout)
    df_ps.wait()
  except CalledProcessError:
    raise Exception("Cannot check {}: No such disk".format(disk_id))
  try:
    output_str = output.decode("utf-8")
  except AttributeError:
    output_str = output
  print(output_str)
  entries = output_str.split()
  # Assumes disk usage percentage is 5th entry of 6
  if len(entries) != 6:
    raise Exception("Unexpected df output. Likely multiple disks identified"
                    " by: {}".format(disk_id))
  usage = entries[4].strip()
  if usage[-1] != "%":
    raise Exception("Unexpected df output. Usage percentage not in 5th col")
  usage = float(usage.strip("%"))
  return usage

def main(sample_id, sample_path, config,  config_path, pipeline):
  disk_id = config["output-disk-identifier"]
  threshold = float(config["free-disk-threshold"])
  disk_usage = _disk_usage(disk_id)
  start_limit = int(config["start-limit"])
  rerun_limit = int(config["rerun-limit"])
  if disk_usage >= threshold:
    print("Disk usage ({}) above threshold of {}".format(disk_usage,
                                                         threshold))
    print("Pausing analysis")
    return
  for index, task in enumerate(pipeline):
    curr_task = task(sample_id, sample_path, config, config_path)
    task_index = str(index).zfill(len(str(len(pipeline))))
    task_complete_flag = "{}/task_{}_finish.flag".format(curr_task.log_dir,
                                                         task_index)
    task_start_flag = "{}/task_{}_begin.flag".format(curr_task.log_dir,
                                                     task_index)
    # if not completed
    if not path.isfile(task_complete_flag):
      # if started before
      if path.isfile(task_start_flag):
        rerun_indices = curr_task.check()
        # if need to rerun
        if rerun_indices:
          print("{} jobs failed".format(len(rerun_indices)))
          with open(task_start_flag, 'r') as start_file:
            times_started = int(start_file.readline().strip())
          if times_started == start_limit:
            reason = "Aborting analysis since current step failed too many times ({}) \n\n{} jobs failed".format(times_started, len(rerun_indices))
            print(reason)
            if config['mail-notify'][0].lower() == 'y':
              msg = ["Aborted", str(curr_task.__class__.__name__), reason, rerun_indices]
              stat = curr_task.mail_notify(msg)
              print(stat)
            break
          elif len(rerun_indices) > rerun_limit:
            reason = "Aborting analysis since too many indices failed ({})".format(len(rerun_indices))
            print(reason)
            if config['mail-notify'][0].lower() == 'y':
              msg = ["Aborted", str(curr_task.__class__.__name__), reason, rerun_indices]
              stat = curr_task.mail_notify(msg)
              print(stat)
            break
          with open(task_start_flag,'w') as start_file:
            start_file.write("{}\n".format(times_started+1))
          curr_task.run(rerun_indices)
          break
        # else job is successful
        else:
          #curr_task.clean_up()
          open(task_complete_flag, 'w').close()
      # else job start for first time
      else:
        with open(task_start_flag, 'w') as start_file:
          start_file.write("1\n")
        curr_task.run()
        break

if __name__ == "__main__":
  pass
