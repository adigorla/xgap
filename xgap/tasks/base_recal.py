#!/ifshome/agorla/data_bucket/apps/python3.7.4/bin/python3

"""Generate Base Recalibration tables"""

try:
  from subprocess import run
except ImportError:
  # Python 2.7 compatibility
  from subprocess import call
  run = call

from os import path, environ, fsync
from sys import stdout, argv
from sys import exit as sysexit
from time import time
from pysam import Samfile, index, AlignmentFile

if 'SLURM_ARRAY_TASK_ID' in  environ:
    taskid = int(environ['SLURM_ARRAY_TASK_ID'])
elif 'SGE_TASK_ID' in environ:
    taskid = int(environ['SGE_TASK_ID'])
elif 'PBS_ARRAYID' in environ:
    taskid = int(environ['PBS_ARRAYID'])

if __name__ == "__main__":
  try:
    INDEX_STR = str(argv[12])
  except IndexError:
    INDEX_STR = str(taskid - 1)
  JAVA_DIR = argv[11]

def base_recalibrator(gatk_jar, in_path, out_path, ref_fa, known_sites,
                      interval_path, n_cthreads, log_output=stdout,
                      bqsr_table=None):
  """GATK BaseRecalibrator
  Generates BQSR Table from input BAM
  Args:
    gatk_jar: Path to GATK jar file
    in_path: Input BAM path
    out_path: Output BQSR Table path
    ref_fa: Path to reference fasta file
    known_sites: List of known_sites paths (such as dbsnp and gold_indel)
    interval_path: Path to interval file for input BAM
    log_output: File handle for log file.
    n_cthreads: Number of computing threads
    bqsr_table: Path to BQSR table if using existing covariates
  """
  cmd = [JAVA_DIR, "-Xmx4g", "-Xms512m","-Djava.awt.headless=true", "-jar", gatk_jar,
         "-T", "BaseRecalibrator",
         "-R", ref_fa,
         "-L", interval_path,
         "-I", in_path,
         "-o", out_path,
         "-nct", n_cthreads]
  for sites in known_sites:
    cmd.append("-knownSites")
    cmd.append(sites)
  # If doing second pass to analyze covariation after recal
  if bqsr_table:
    cmd.append("-BQSR")
    cmd.append(bqsr_table)
  start = time()
  run(cmd, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("BaseRecalibrator completed in {} seconds\n".format(end-start))
  log_output.flush()
  fsync(log_output.fileno())

def _cmp_lt(alignment_1, alignment_2):
  """Equivalent to alignment_1 < alignment_2, using reference position"""
  if alignment_1.reference_id == alignment_2.reference_id:
    return alignment_1.reference_start < alignment_2.reference_start
  else:
    return alignment_1.reference_id < alignment_2.reference_id

def merge_two_bams(input_path_1, input_path_2, output_path):
  """Merges two sorted BAMs"""
  with Samfile(input_path_1, 'rb') as input_file_1,\
       Samfile(input_path_2, 'rb') as input_file_2,\
       Samfile(output_path, 'wb', header = dict(input_file_1.header)) as output_file:
    read_1 = next(input_file_1, None)
    read_2 = next(input_file_2, None)
    while read_1 and read_2:
      if _cmp_lt(read_1, read_2):
        output_file.write(read_1)
        read_1 = next(input_file_1, None)
      else:
        output_file.write(read_2)
        read_2 = next(input_file_2, None)
    while read_1:
      output_file.write(read_1)
      read_1 = next(input_file_1, None)
    while read_2:
      output_file.write(read_2)
      read_2 = next(input_file_2, None)


def main(gatk_jar, sample_id, out_dir, ref_fa, known_sites_str, interval_dir,
         n_cthreads, n_regions, log_prefix, sambamba):
  log_path = "{}.{}.log".format(log_prefix, INDEX_STR)
  log_output = open(log_path, 'w')
  region_num = INDEX_STR.zfill(len(n_regions))
  region_name = "region_{}".format(region_num)
  region_bam = "{}/Dedup/{}_{}.dedup.bam".format(out_dir, sample_id,
                                                      region_name)
  chri_bam = "{}/Dedup/{}_chri.dedup_{}.bam".format(out_dir, sample_id,
                                                    region_name)
  known_sites = [site.strip() for site in known_sites_str.split(',')]
  start = time()
  input_bam = "{}/Recal/{}_{}.dedup.bam".format(out_dir, sample_id, region_name)
  
  '''
  #merge and index bam files
  cmd1 = [sambamba,"merge",
         "--nthreads=1",
         "--compression-level=3"]
  cmd1 += [input_bam]
  cmd1 += [region_bam, chri_bam]
  log_output.write("Merging chrI reads\n")
  run(cmd1, stdout=log_output, stderr=log_output)
  log_output.write("Merged chrI reads\n")
  '''
  

#merge chri                                                                                                                                                                       
  input_bam = "{}/Recal/{}_{}.dedup.bam".format(out_dir, sample_id, region_name)
  log_output.write("Merging chrI reads\n")
  log_output.flush()
  fsync(log_output.fileno())
  merge_two_bams(region_bam, chri_bam, input_bam)
  log_output.write("Merged chrI reads\n")
  log_output.flush()
  fsync(log_output.fileno())  
#index BAM files                                                                                                                                                                  
  try:
    AlignmentFile(input_bam, "rb").check_index()
  except:
    log_output.write("Generating index file")
    log_output.flush()
    fsync(log_output.fileno())
    index(input_bam,"{}.bai".format(input_bam))

  # Generate BQSR table
  log_output.write("Generating BQSR table\n")
  log_output.flush()
  fsync(log_output.fileno())
  interval_path = "{}/{}.m2bp.intervals".format(interval_dir, region_name)
  bqsr_path = "{}/Recal/{}_{}.dedup.bqsr.csv".format(out_dir, sample_id,
                                                     region_num)
  base_recalibrator(gatk_jar, input_bam, bqsr_path, ref_fa, known_sites,
                    interval_path, n_cthreads, log_output)
  end = time()
  log_output.write("Merge chrI and BQSR completed in "
                   "{} seconds\n".format(round(end-start, 2)))
  log_output.close()
  sysexit()

if __name__ == "__main__":
  gatk_jar = argv[1]
  sample_id = argv[2]
  out_dir = argv[3]
  ref_fa = argv[4]
  known_sites_str = argv[5]
  interval_dir = argv[6]
  n_cthreads = argv[7]
  n_regions = argv[8]
  log_prefix = argv[9]
  sambamba  = argv[10]
  main(gatk_jar, sample_id, out_dir, ref_fa, known_sites_str, interval_dir,
       n_cthreads, n_regions, log_prefix, sambamba)
