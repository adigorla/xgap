#!/ifshome/agorla/data_bucket/apps/python3.7.4/bin/python3

"""Merges BQSR tables for each region"""
"""IMP NOTE: This BQSR report merge script only works with GATK 3.7 or above"""
try:
  from subprocess import run
except ImportError:
  # Python 2.7 compatibility
  from subprocess import call
  run = call

from os import fsync
from sys import stdout, argv
from sys import exit as sysexit
from time import time

def bqsr_gather(gatk_jar, in_paths, out_path, log_output=stdout):
  """GATK GatherBqsrReports
  Combines BQSR Tables
  Args:
    gatk_jar: Path to GATK jar file
    in_paths: List of paths to each BQSR Table to be merged
    out_path: Path to output merged BQSR Table
    log_output: File handle for log file
  """
  # Tool name subject to change.
  output=checkout([JAVA_DIR, "-Xmx4g", "-Xms512m","-Djava.awt.headless=true", "-jar", gatk_jar, "--version"])
  charstr=output.decode('ASCII')
  gatk-ver="4"
  for i, c in enumerate(charstr):
    if c.isdigit():
        gatk-ver=(charstr[i])
        break
  if gatk-ver="3":
    gather_tool = "org.broadinstitute.gatk.tools.GatherBqsrReports"
  elif gatk-ver="4":
    gather_tool = "GatherBQSRReports"
  cmd = [JAVA_DIR, "-Xmx4g", "-Djava.awt.headless=true", "-cp", gatk_jar, gather_tool]
  for in_path in in_paths:
    cmd.append("I={}".format(in_path))
  cmd.append("O={}".format(out_path))
  start = time()
  run(cmd, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("GatherBQSR completed in {} seconds\n".format(round(end-start, 2)))
  log_output.flush()
  fsync(log_output.fileno())

def main(gatk_jar, sample_id, out_dir, n_regions, log_path):
  in_paths = []
  for index in range(int(n_regions)):
    region_name = str(index).zfill(len(n_regions))
    in_bqsr = "{}/Recal/{}_{}.dedup.bqsr.csv".format(out_dir, sample_id,
                                                     region_name)
    in_paths.append(in_bqsr)
  out_path = "{}/Recal/{}.dedup.bqsr.csv".format(out_dir, sample_id)
  log_output = open(log_path, 'w')
  log_output.write("Merging BQSR tables\n")
  log_output.flush()
  fsync(log_output.fileno())
  bqsr_gather(gatk_jar, in_paths, out_path, log_output)
  log_output.close()
  sysexit()

if __name__ == "__main__":
  gatk_jar = argv[1]
  sample_id = argv[2]
  out_dir = argv[3]
  n_regions = argv[4]
  log_path = argv[5]
  JAVA_DIR = argv[6]
  main(gatk_jar, sample_id, out_dir, n_regions, log_path)
