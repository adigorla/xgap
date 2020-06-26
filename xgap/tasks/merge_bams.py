#!/ifshome/agorla/data_bucket/apps/python3.7.4/bin/python3

"""Methods for merging sorted BAM files and dedupping"""

try:
  from subprocess import run
except ImportError:
  # Python 2 compatibility
  from subprocess import call
  run = call

from glob import glob
from os import remove, path, fsync
from sys import stdout, argv
from sys import exit as sysexit
from time import time

from pysam import Samfile, merge

def _cmp_lt(alignment_1, alignment_2):
  """Equivalent to alignment_1 < alignment_2, using reference position"""
  if alignment_1.reference_id == alignment_2.reference_id:
    return alignment_1.reference_start < alignment_2.reference_start
  else:
    return alignment_1.reference_id < alignment_2.reference_id

def _load_region(index, region_dir, n_regions, seq_dict):
  """Load region file into dict structure"""
  region = {}
  index_str = str(index).zfill(len(str(n_regions)))
  region_path = "{}/region_{}.m0bp.intervals".format(region_dir, index_str)
  with open(region_path, 'r') as region_file:
    for line in region_file:
      reference_id, interval = line.strip().split(":")
      lower_bound, upper_bound = [int(i) for i in interval.split("-")]
      region[seq_dict[reference_id]] = [lower_bound, upper_bound]
  return region

def merge_interregional_bams(input_paths, output_path, regions_dir):
  """Merges regional sorted BAMs
  Should be used for regional BAMs that are sorted and have defined intervals.
  Written to work with XGAP pipeline only, specifically using config to
  find interval files.
  Args:
    input_paths: List of regional BAMs, in proper order
    output_path: Path to output merged BAM
    regions_dir: Path to directory with interval files
    config: Config dict for XGAP
  """
  n_regions = len(input_paths)
  sample = Samfile(input_paths[0], 'rb')
  sample_header = dict(sample.header)
  sample.close()
  seq_dict = {}
  for i, seq in enumerate(sample_header["SQ"]):
    seq_dict[seq["SN"]] = i
  #regions_dir = "{}/regions/{}/".format(config['output-dir'], n_regions)
  with Samfile(output_path, 'wb', header=sample_header) as out_file:
    for index, input_path in enumerate(input_paths):
      # Object Info: region['ref_id'] = [lower_bound, upper_bound]
      region = _load_region(index, regions_dir, n_regions, seq_dict)
      started = False
      with Samfile(input_path, 'rb') as in_file:
        for alignment in in_file:
          pos = alignment.reference_start
          ref_id = alignment.reference_id
          if region[ref_id][0] <= pos <= region[ref_id][1]:
            if not started:
              started = True
            out_file.write(alignment)
          else:
            # Reached end of interval, move to next file
            if started:
              break

def main(sample_id, out_dir, n_regions, regions_dir, log_path):
  log_output = open(log_path , 'w')
  complete_start = time()

  # Merge regional BAMs
  log_output.write("Merging regional BAMs\n")
  log_output.flush()
  fsync(log_output.fileno())
  in_bams = []
  for index in range(int(n_regions)):
    region_name = str(index).zfill(len(n_regions))
    in_bam = "{}/bam/{}_{}.dedup.recal.bam".format(out_dir, sample_id,
                                                   region_name)
    in_bams.append(in_bam)
  out_bam = "{}/bam/{}.dedup.recal.bam".format(out_dir, sample_id)
  start = time()
  merge_interregional_bams(in_bams, out_bam, regions_dir)
  end = time()
  if path.isfile(out_bam):
    log_output.write("Merged regional BAMs in {} seconds\n".format(end - start))
    log_output.flush()
    fsync(log_output.fileno())

#TODO: -use sambamba to merge bams
  # Merge qcfail BAMs
  log_output.write("Merging qcfail BAMs\n")
  qcfail_dir = "{}/bam/qcfail".format(out_dir)
  in_bams = glob("{}/*.bam".format(qcfail_dir))
  out_bam = "{}/../{}_qcfail.bam".format(qcfail_dir, sample_id)
  bam_list_path = "{}/bam_list.txt".format(qcfail_dir)
  with open(bam_list_path, "w") as bam_list:
    for bam in in_bams:
      bam_list.write("{}\n".format(bam))
  merge('-f', "-b", bam_list_path, out_bam)
  if path.isfile(out_bam):
    log_output.write("Merged qcfail BAMs\n")
    log_output.flush()
    fsync(log_output.fileno())
    remove(bam_list_path)

  # Merge unmapped BAMs
  log_output.write("Merging unmapped BAMs\n")
  log_output.flush()
  fsync(log_output.fileno())
  unmapped_dir = "{}/bam/unmapped".format(out_dir)
  in_bams = glob("{}/*.bam".format(unmapped_dir))
  out_bam = "{}/../{}_unmapped.bam".format(unmapped_dir, sample_id)
  bam_list_path = "{}/bam_list.txt".format(unmapped_dir)
  with open(bam_list_path, "w") as bam_list:
    for bam in in_bams:
      bam_list.write("{}\n".format(bam))
  merge('-f', "-b", bam_list_path, out_bam)
  if path.isfile(out_bam):
    log_output.write("Merged unmapped BAMs\n")
    log_output.flush()
    fsync(log_output.fileno())
    remove(bam_list_path)
  complete_end = time()  
  log_output.write("Merged all BAMs in "
                   "{} seconds\n".format(round(complete_end-complete_start)))
  log_output.close()
  sysexit()

if __name__ == "__main__":
  sample_id = argv[1]
  out_dir = argv[2]
  n_regions = argv[3]
  regions_dir = argv[4]
  log_path = argv[5]
  main(sample_id, out_dir, n_regions, regions_dir, log_path)
