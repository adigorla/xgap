#!/u/local/apps/python/3.7.2/bin/python3

"""Aligns reads and sorts into proper regions and coordinate order"""

from collections import namedtuple
from glob import glob
from os import path, environ, fsync
from subprocess import Popen, PIPE
from sys import stdout, argv
from sys import exit as sysexit
from time import time
from pysam import Samfile, AlignedRead, qualitystring_to_array

if 'SLURM_ARRAY_TASK_ID' in  environ:
    taskid = int(environ['SLURM_ARRAY_TASK_ID'])
elif 'SGE_TASK_ID' in environ:
    if environ['SGE_TASK_ID'] != 'undefined':
        taskid = int(environ['SGE_TASK_ID'])
elif 'PBS_ARRAYID' in environ:
    taskid = int(environ['PBS_ARRAYID'])

__all__ = ['align_and_sort']

# BUILT FOR SGE 
#   To change for other job managers, make sure this variable indicates which FASTQ chunk
#   this job will be handling ( should be int in [0, num_chunks) ) 
if __name__ == "__main__":
  try:
    index_str = str(argv[9])
  except IndexError:
    index_str = str(taskid - 1)

# Python 2/3 compatibility
try:
  range = xrange
except NameError:
  pass

# BWA algorithm to use. Currently written for bwa mem
_ALGORITHM = "mem"

# Region object for regions dict. Described in _load_regions()
_Region = namedtuple("Region", ["name", "lower_bound", "upper_bound"])

#CIGAR values
_CIGAR_OPERATIONS = {}
for _val, _char in enumerate("MIDNSHP=XB"):
  _CIGAR_OPERATIONS[_char] = _val

def _gen_bwa_cmd(bin_bwa, n_threads, ref_fa, in_fq, sample_id, platform_id):
  """Generate bwa command with given arguments
  "-M" option added to command for Picard compatibility.
  Args:
    bin_bwa: Full path to bwa executable
    n_threads: Number of threads to use for alignment process
    ref_fa: Full path to reference fasta file
    in_fq: List of input fastq file(s).
    sample_id: Unique identifier for this sample (for output dir)
    platform_id: sequencing platform info
  Returns:
    bwa_cmd: String that can be run in command line.
      (Ex. "bwa mem -M -t 1 ref.fa test.fq")
  """
  read_group = "\"@RG\\tID:{sid}\\tPL:{pid}\\tSM:{sid}\"".format(sid=sample_id, pid=platform_id)
  bwa_cmd = [bin_bwa, _ALGORITHM, "-M",
             "-t", n_threads,
             "-R", read_group]
  if len(in_fq) == 1:
    bwa_cmd.append("-p")
  bwa_cmd.append(ref_fa)
  for fastq in in_fq:
    bwa_cmd.append(fastq)
  return " ".join(bwa_cmd)

def _load_regions(interval_paths):
  """Loads region information into a convenient dict structure
  Args:
    interval_paths: List of sorted interval file paths to load
  Returns:
    regions: Dict storing region information.
      regions['seq'] = ["Region" objects corresponding to regions in 'seq']
      "Region" object is a named tuple with members:
        index: Index of corresponding region
        lower_bound: Lower bound of region's interval in 'seq'
        upper_bound: Upper bound of region's interval in 'seq'
      Moreover, this list of Region objects is sorted to allow binary search
  """
  regions = {}
  for region_path in interval_paths:
    region_name = path.splitext(path.splitext(path.basename(region_path))[0])[0] #region_XXXX
    with open(region_path, 'r') as region_file:
      for line in region_file:
        chromosome, bounds = line.strip().split(":")
        lower_bound, upper_bound = [int(bound) for bound in bounds.split("-")]
        region = _Region(region_name, lower_bound, upper_bound)
        if chromosome in regions:
          regions[chromosome].append(region)
        else:
          regions[chromosome] = [region]
  return regions

def _get_regions(read, regions, seq_dict, log_output):
  """Return region name corresponding to given chromosome and position
  Args:
    read: AlignedRead()
    regions: Dict structure produced by _load_regions()
    seq_dict: Dictionary mapping reference_id to string representation
  Returns:
    region_names: list of strings corresponding to regions encompassing position
  """
  chromosome = seq_dict[read.reference_id]
  position = read.reference_start + 1
  region_names = []
  curr_regions = regions[chromosome]
  first = 0
  last = len(curr_regions) - 1
  while first <= last:
    midpoint = (first + last) // 2
    curr_region = curr_regions[midpoint]
    if curr_region.lower_bound <= position <= curr_region.upper_bound:
      region_names.append(curr_region.name)
      index = midpoint + 1
      while index < len(curr_regions):
        curr_region = curr_regions[index]
        if curr_region.lower_bound <= position <= curr_region.upper_bound:
          region_names.append(curr_region.name)
          index += 1
        break
      index = midpoint - 1
      while index >= 0:
        curr_region = curr_regions[index]
        if curr_region.lower_bound <= position <= curr_region.upper_bound:
          region_names.append(curr_region.name)
          index -= 1
        break
      return sorted(region_names)
    else:
      if position < curr_region.lower_bound:
        last = midpoint - 1
      else:
        first = midpoint + 1
  err_msg = ("Read at ({}, {}) does not fall into "
             "any region".format(chromosome, position))
  log_output.write("ERROR: {}\n".format(err_msg))
  raise Exception(err_msg)

def _sort_regional_reads(regional_reads, regions, output_dir, basename, header,
                         log_output):
  """Sorts reads in dict struct by coordinate and writes to BAM files
  Uses counting sort with the intuition that most positions in a range will
  map to a non-empty bucket, where a bucket holds reads mapped to a certain
  position.
  Args:
    regional_reads: Complicated dict struct
      ex. regional_reads[region][chromosome][position] returns list of reads
          at this position in this chromosome in this region
    regions: Regions dict used to sort reads during alignment. Produced by
      _load_regions()
    output_dir: self-explanatory
    basename: Prefix for output BAM files
    header: Dict representation of SAM header (described in pysam docs)
    log_output: File object for log.
  Returns:
    output_paths: List of output file paths
  """
  output_paths = []
  region_intervals = {}
  sequence_order = [chrom["SN"] for chrom in header["SQ"]]
  region_names = set(["chri"])
  # Get region intervals
  for seq_id in sequence_order:
    for interval in regions[seq_id]:
      region_names.add(interval.name)
      if interval.name in region_intervals:
        region_intervals[interval.name].append([seq_id, interval.lower_bound,
                                                interval.upper_bound])
      else:
        region_intervals[interval.name] = [[seq_id, interval.lower_bound,
                                            interval.upper_bound]]
  #for region in regional_reads.keys():
  for region in region_names:
    output_path = "{}/{}/{}.bam".format(output_dir, region, basename)
    output_file = Samfile(output_path, 'wb', header=header)
    if region == "chri":
      for chromosome in sequence_order:
        if chromosome in regional_reads[region]:
          for position in sorted(regional_reads[region][chromosome].keys()):
            for read in regional_reads[region][chromosome].pop(position):
              output_file.write(read)
    elif region in regional_reads:
      for entries in region_intervals[region]:
        chromosome = entries[0]
        low, high = entries[1], entries[2]
        if chromosome in regional_reads[region]:
          for position in range(low, high + 1):
            if position in regional_reads[region][chromosome]:
              for read in regional_reads[region][chromosome].pop(position):
                output_file.write(read)
        else:
          log_output.write("No reads mapping to "
                           "{} {}\n".format(region, chromosome))
          log_output.flush()
          fsync(log_output.fileno())
    output_paths.append(output_path)
    output_file.close()
  return output_paths

def _add_read_to_dict(read, region, regional_reads, seq_dict):
  """Adds alignment to hash structure for sorting
  Args:
    read: AlignedRead
    region: Region to add read to
    regional_reads is a dict storing reads by region.
    seq_dict: Dict mapping string for ref ID to index.
  """
  position = read.reference_start + 1
  chromosome = seq_dict[read.reference_id]
  if region in regional_reads:
    if chromosome in regional_reads[region]:
      if position in regional_reads[region][chromosome]:
        regional_reads[region][chromosome][position].append(read)
      else:
        regional_reads[region][chromosome][position] = [read]
    else:
      regional_reads[region][chromosome] = {position: [read]}
  else:
    regional_reads[region] = {chromosome: {position: [read]}}

def _add_line_to_header(line, header):
  """Parse header lines and adds to dict representation of header"""
  entries = line.strip().split('\t')
  first_key = entries[0].strip('@')
  second_dict = {}
  for entry in entries[1::]:
    second_key = entry[0:2]
    second_val = entry[3::]
    second_dict[second_key] = second_val
  if first_key in header:
    if first_key == "HD":
      for key in second_dict:
        header[first_key][key] = second_dict[key]
    else:
      header[first_key].append(second_dict)
  else:
    if first_key == "HD":
      header[first_key] = second_dict
    else:
      header[first_key] = [second_dict]

def _string_to_aligned_segment(line, seq_dict, log_output):
  """Converts SAM record in string format to pysam AlignedRead
  Args:
    line: String of SAM record
    seq_dict: Dictionary mapping reference ID to reference ID index
    log_output: Handle for outputting log information
  Returns:
    aligned_segment: pysam AlignedRead class with values from 'line'
  """
  line = line.strip().split()
  aligned_segment = AlignedRead()
  aligned_segment.query_name = line[0]
  aligned_segment.flag = int(line[1])
  if line[2] != "*":
    aligned_segment.reference_id = seq_dict[line[2]]
    aligned_segment.reference_start = int(line[3]) - 1
    aligned_segment.mapping_quality = int(line[4])
  cigartuples = []
  pos = ""
  for symbol in line[5]:
    if symbol.isdigit():
      pos += symbol
    elif symbol == "*":
      continue
    else:
      cigartuples.append((_CIGAR_OPERATIONS[symbol], int(pos)))
      pos = ""
  aligned_segment.cigartuples = cigartuples
  if line[6] == "=":
    aligned_segment.next_reference_id = seq_dict[line[2]]
  elif line[6] != "*":
    aligned_segment.next_reference_id = seq_dict[line[6]]
  aligned_segment.next_reference_start = int(line[7]) - 1
  aligned_segment.template_length = int(line[8])
  aligned_segment.query_sequence = line[9]
  aligned_segment.query_qualities = qualitystring_to_array(line[10])
  for field in line[11::]:
    tag, tag_type, val = field.split(":")
    if tag_type == "i":
      val = int(val)
    elif tag_type == "f":
      val = float(val)
    elif tag_type == "H":
      val = bytearray.fromhex(val)
    elif tag_type == "B":
      val = [int(i) for i in val.split(",")]
    elif not (tag_type == "A" or tag_type == "Z"):
      err_msg = "Optional Ttag type '{}' not recognized".format(tag_type)
      log_output.write("ERROR: {}\n".format(err_msg))
      raise Exception(err_msg)
    aligned_segment.set_tag(tag, val, value_type=tag_type)
  return aligned_segment

def _read_sorter(rlist, unmapped_file, qcfail_file, regions, regional_reads, seq_dict, log_output):
  """                                                                                                                                                                                
  Function to sort reads with the same read name                                                                                                                                     
  """
  prime_reads = []
  r1, r1r, r2, r2r = None, None, None, None

  for uaread in rlist:
    if uaread.is_unmapped:
      unmapped_file.write(uaread)
      continue
    if uaread.is_qcfail:
      qcfail_file.write(uaread)
      continue
    reg = _get_regions(uaread, regions, seq_dict, log_output)
    #condtion to find primary reads                                                                                                                                                  
    #if not primary assign read to bwa assigned pos/region
    if uaread.is_secondary or uaread.is_supplementary:
      # if non-primary read is mapped to regional boundary assign to chrI                                                                                                            
      if len(reg) > 1:
        _add_read_to_dict(uaread, "chri", regional_reads, seq_dict)
      else:
        _add_read_to_dict(uaread, reg[0], regional_reads, seq_dict)
    #if primary read save to new list for further sorting                                                                                                                            
    else:
      prime_reads.append(uaread)
  
  #check if a PAIR of primary reads is identified                                                                                                                                    
  if len(prime_reads) == 2:
    r1 = prime_reads[0]
    r2 = prime_reads[1]
    r1r = _get_regions(r1, regions, seq_dict, log_output)
    r2r = _get_regions(r2, regions, seq_dict, log_output)
    #checks if both primary reads are assigned to same region + chr and NOT in regional bondary                                                                                      
    if len(r1r) == 1 and len(r2r) == 1 and seq_dict[r1.reference_id] == seq_dict[r2.reference_id] and r1r == r2r:
      _add_read_to_dict(r1, r1r[0], regional_reads, seq_dict)
      _add_read_to_dict(r2, r2r[0], regional_reads, seq_dict)
    else:
      _add_read_to_dict(r1, "chri", regional_reads, seq_dict)
      _add_read_to_dict(r2, "chri", regional_reads, seq_dict)
  else:
    for pread in prime_reads:
      reg = _get_regions(pread, regions, seq_dict, log_output)
      if len(reg) > 1:
        _add_read_to_dict(pread, "chri", regional_reads, seq_dict)
      else:
        _add_read_to_dict(pread, reg[0], regional_reads, seq_dict)

def _sort_alignments_by_region(bwa_cmd, regions, output_dir, basename,
                               log_output=stdout):
  """Sort bwa alignments into specified regions
  Sorts alignments from bwa command into separate files for each specified
  region. Sorts unmapped and qc fail reads into separate files. Paired reads
  that fall into different regions are added to chrI file.
  Args:
    bwa_cmd: String for full bwa command (Ex. "bwa mem -M -t 1 ref.fa test.fq")
    regions: A dict storing region information.
      regions['seq'] = [list of regions that encompass a part of 'seq']
    output_dir: Directory storing region subdirectories that store alignments
    basename: File prefix for output BAM files.
    log_output: Handle for log output
  Returns:
    output_paths: List of all output files.
  """
  regional_reads, header, seq_dict = {}, {}, {}
  curr_reads = []
  read = None
  unmapped_path = "{}/unmapped/{}.bam".format(output_dir, basename)
  qcfail_path = "{}/qcfail/{}.bam".format(output_dir, basename)
  log_output.write("executing BWA\n")
  start = time()
  bwa_process = Popen(bwa_cmd, stdout=PIPE, stderr=log_output, shell=True)

  for line in bwa_process.stdout:
    line = line.decode("utf-8")
    if line[0] == '@':
      _add_line_to_header(line, header)
      continue
    if not seq_dict:
      for index, seq in enumerate(header["SQ"]):
        header["SQ"][index]["LN"] = int(header["SQ"][index]["LN"])
        seq_dict[seq["SN"]] = index
        seq_dict[index] = seq["SN"]
      unmapped_file = Samfile(unmapped_path, 'wb', header=header)
      qcfail_file = Samfile(qcfail_path, 'wb', header=header)
    read = _string_to_aligned_segment(line, seq_dict, log_output)
    #conditionals to group bwa output by read names and sort one group at a time                                                                                                     
    if not curr_reads:
      curr_reads.append(read)
      continue
    if curr_reads[0].query_name == read.query_name:
      curr_reads.append(read)
      continue
    else:
      _read_sorter(curr_reads, unmapped_file, qcfail_file, regions, regional_reads, seq_dict, log_output)
      curr_reads.clear()
      curr_reads.append(read)
  if curr_reads:
    _read_sorter(curr_reads, unmapped_file, qcfail_file, regions, regional_reads, seq_dict, log_output)
    curr_reads.clear()

  _add_line_to_header("@HD\tSO:coordinate", header)
  output_paths = _sort_regional_reads(regional_reads, regions, output_dir, basename,
                                      header, log_output)
  end = time()
  log_output.write("Sorted {} regions in {} "
                   "seconds\nDone\n".format(len(output_paths), (end - start)))
  log_output.flush()
  fsync(log_output.fileno())
  output_paths.append(unmapped_path)
  output_paths.append(qcfail_path)
  unmapped_file.close()
  qcfail_file.close()
  return output_paths

def _get_interval_paths(regions_dir, extension):
  """Returns paths of interval files in regions_dir with given
     extension in a sorted list"""
  return sorted(glob("{}/*{}".format(regions_dir, extension)))

def align_and_sort(bin_bwa, n_threads, ref_fa, platform_id, in_fq, sample_id, output_dir,
                   interval_paths, log_output=stdout):
  """Aligns in_fq using bwa and outputs sorted regional BAM chunks to output dir
  Args:
    bin_bwa: Path to bwa executable
    n_threads: Available threads to parallelize bwa
    ref_fa: Path to reference FASTA
    in_fq: List of input fastq files.
    sample_id: Unique identifier for sample.
    output_dir: Directory to place sorted files
    interval_paths: List of sorted interval file paths defining regions
    log_output: Handle to write log output to.
    platform_id: sequencing platform info
  Returns:
    output_paths: List of output BAM file paths.
  """
  bwa_cmd = _gen_bwa_cmd(bin_bwa, n_threads, ref_fa, in_fq, sample_id, platform_id)
  regions = _load_regions(interval_paths)
  basename = path.basename(in_fq[0]).split('.')[0]
  log_output.write("Strating alignment_by_region\n")
  return _sort_alignments_by_region(bwa_cmd, regions, output_dir, basename,
                                    log_output)

def split_bam_into_regions(in_path, out_prefix, interval_dir, log_output):
  """Splits a sorted BAM file into regions"""
  regions = _load_regions(_get_interval_paths(interval_dir, "m3bp.intervals"))
  region_names = set()
  for seq_id in regions.keys():
    for interval in regions[seq_id]:
      region_names.add(interval.name)
  region_names = sorted(region_names)
  open_bams = {}
  bams_opened = set()
  with Samfile(in_path, 'rb') as in_bam:
    seq_dict = {}
    for index, seq in enumerate(in_bam.header["SQ"]):
      seq_dict[seq["SN"]] = index
      seq_dict[index] = seq["SN"]
    for read in in_bam:
      regions_written = set()
      for region in _get_regions(read, regions, seq_dict, log_output):
        regions_written.add(region)
        if region not in open_bams:
          open_bams[region] = Samfile("{}_{}.bam".format(out_prefix, region),
                                      'wb', template=in_bam)
          if region in bams_opened:
            err_msg = "Input BAM not sorted"
            log_output.write("ERROR: {}\n".format(err_msg))
            raise Exception(err_msg)
          bams_opened.add(region)
        open_bams[region].write(read)
      for region in list(open_bams):
        if region not in regions_written:
          open_bams.pop(region).close()
    for region in region_names:
      if region not in bams_opened:
        Samfile("{}_{}.bam".format(out_prefix, region), 'wb', template=in_bam).close()

def main(bin_bwa, n_threads, ref_fa, sample_id, out_dir, interval_dir, log_prefix, platform_id):
#TODO: check if sambamba will be better/faster at sorting!
  log_path = "{}.{}.log".format(log_prefix, index_str)
  log_output = open(log_path, 'w')
  fastq_prefix = "{}/Fastq/{}".format(out_dir, sample_id)
  in_fq = ["{}_{}_{}.fq.gz".format(fastq_prefix, index_str, "R1"),
           "{}_{}_{}.fq.gz".format(fastq_prefix, index_str, "R2")]
  output_dir = "{}/Bwa/".format(out_dir)
  interval_paths = _get_interval_paths(interval_dir,"m3bp.intervals")
  log_output.write("Strating align_and_sort\n")
  start = time()
  output_paths = align_and_sort(bin_bwa, n_threads, ref_fa, platform_id, in_fq, sample_id,
                                output_dir, interval_paths, log_output)
  end = time()
  if output_paths:
    log_output.write("Alignment and sorting completed in {} seconds\n".format(round(end-start, 2)))
  else:
    log_output.write("ERROR: No output files generated\n")
  log_output.close()
  sysexit()

if __name__ == "__main__":
  # path to bwa binary
  bin_bwa = argv[1]
  # number of available threads
  n_threads = argv[2]
  # path to reference fasta
  ref_fa = argv[3]
  sample_id = argv[4]
  # sample output dir
  out_dir = argv[5]
  # directory with interval files
  interval_dir = argv[6]
  # log_prefix
  log_prefix = argv[7]
  #sequencing platform
  platform_id = argv[8]
  main(bin_bwa, n_threads, ref_fa, sample_id, out_dir, interval_dir, log_prefix, platform_id)
