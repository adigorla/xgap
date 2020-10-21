#!/u/local/apps/python/3.7.2/bin/python3

"""Converts input sequence files to split and compressed FASTQ files"""

from glob import glob
from gzip import open as gzopen
from sys import stdin, stdout, argv
import sys
from pysam import Samfile
from time import time

__all__ = ['split_input', 'split_sam_stream', 'split_alignment_file',
           'split_paired_fastq', 'split_interleaved_fastq']

#Global constants for development
_RECORDS_PER_FILE = 750000 # ~60MB
_UPDATE_COUNT = 10000000 # Records processed to update log file
_COMPRESS_LVL = 5
_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

# Filter options
_FILTER_QC = True
_FILTER_NON_PRIMARY = True
_FILTER_UNPAIRED = True

# SAM flag values
_IS_PAIRED = 1
_IS_REVERSE = 16
_IS_READ1 = 64
_IS_READ2 = 128
_IS_NOT_PRIMARY_ALIGNMENT = 256
_IS_QC_FAIL = 512

class _RecordWriter(object):
  """Interface for writing records to split paired-end FASTQ files
    Formats records to FASTQ format and writes to files in chunks of
    roughly 60 MB for parallelization of downstream analysis. Unpaired
    records are written to one file if not filtered.
  Attributes:
    fastq_prefix: Full path prefix for FASTQ files generated. Index string
      and file extension will be appended to this.
    fastq_index: Current index of chunk being written to.
    fastq_records: Number of records written to current chunk.
    fastq_file_1: Current chunk for first-in-pair records.
    fastq_file_2: Current chunk for second-in-pair records.
    unpaired_file: File to write all unpaired records to.
  """
  def __init__(self, fastq_prefix, log_output=stdout, output_queue=None):
    """Inits _RecordWriter for given FASTQ prefix
    Args:
      output_queue: Queue() object from queue builtin. If provided, output
        chunks will be added as they are created. Multithreading compatible.
    """
    self.fastq_prefix = fastq_prefix
    self.fastq_index = -1
    self.fastq_records = 0
    self.fastq_file_1 = None
    self.fastq_file_2 = None
    self.unpaired_file = None
    self.output_paths = []
    self.log_output = log_output
    self.output_queue = output_queue
    self.update_fastq_index()

  def update_fastq_index(self):
    """Closes current chunks for FASTQ files and creates next chunks"""
    if self.fastq_file_1:
      self.fastq_file_1.close()
      self.fastq_file_2.close()
    self.fastq_index += 1
    self.fastq_records = 0
    #index_str = str(self.fastq_index).zfill(5)
    # files won't be sorted but easier for indexed job management
    index_str = str(self.fastq_index)
    fq_path_1 = '{}_{}_{}.fq.gz'.format(self.fastq_prefix, index_str, "R1")
    fq_path_2 = '{}_{}_{}.fq.gz'.format(self.fastq_prefix, index_str, "R2")
    self.fastq_file_1 = gzopen(fq_path_1, 'wt', compresslevel=_COMPRESS_LVL)
    self.fastq_file_2 = gzopen(fq_path_2, 'wt', compresslevel=_COMPRESS_LVL)
    #print("------------- Record writer atributes {} ---------------------------------".format(self.__dict__))
    self.output_paths.append([self.fastq_file_1.name,
                              self.fastq_file_2.name])
    if self.output_queue:
      self.output_queue.put([self.fastq_file_1.name,
                              self.fastq_file_2.name])

  def write_paired_records(self, record_a, record_b, pair_id=True):
    """Writes given records to appropriate FASTQ chunk
    Determines which read is 1st in pair, and whether each sequence maps to
    reverse strand. Writes appropriately formatted records to respective chunk
    Args:
      record_a, record_b: Paired records (pysam.Samfile.fetch() item)
      in no particular order
    """
    if record_a.is_read1 and record_b.is_read2:
      record_1 = record_a
      record_2 = record_b
    elif record_a.is_read2 and record_b.is_read1:
      record_1 = record_b
      record_2 = record_a
    else:
      self.log_output.write("Not proper pair")
      return
    if pair_id:
      qname_1 = '@{}/1'.format(record_1.qname)
      qname_2 = '@{}/2'.format(record_2.qname)
    else:
      qname_1 = '@{}'.format(record_1.qname)
      qname_2 = '@{}'.format(record_2.qname)
    if record_1.is_reverse:
      self.fastq_file_1.write("\n".join([qname_1,
                                         _reverse_complement(record_1.seq),
                                         "+", record_1.qual, ""]))
    else:
      self.fastq_file_1.write("\n".join([qname_1,
                                         record_1.seq, "+",
                                         record_1.qual, ""]))
    if record_2.is_reverse:
      self.fastq_file_2.write("\n".join([qname_2,
                                         _reverse_complement(record_2.seq),
                                         "+", record_2.qual, ""]))
    else:
      self.fastq_file_2.write("\n".join([qname_2,
                                         record_2.seq, "+",
                                         record_2.qual, ""]))
    self.fastq_records += 1
    if self.fastq_records == _RECORDS_PER_FILE:
      self.update_fastq_index()

  def write_unpaired_record(self, record):
    """Writes given record to unpaired FASTQ file and formats appropriately"""
    if not self.unpaired_file:
      unpaired_path = ''.join([self.fastq_prefix, '_unpaired.fq.gz'])
      self.unpaired_file = gzopen(unpaired_path, 'wt',
                                  compresslevel=_COMPRESS_LVL)
      self.output_paths.append(unpaired_path)
    qname = '@{}'.format(record.qname)
    if record.is_reverse:
      self.unpaired_file.write("\n".join(
        [qname, _reverse_complement(record.seq), "+", record.qual, ""]))
    else:
      self.unpaired_file.write("\n".join(
        [qname, record.seq, "+", record.qual, ""]))

class _SimpleRecord(object):
  """Object that stores the record entries needed for a FASTQ record
  Attributes:
    qname: Name of record
    seq: Sequence string
    qual: Quality string
  """
  def __init__(self, qname=None, seq=None, qual=None,
               is_read1=False, is_read2=False):
    """Inits object with FASTQ entries"""
    self.qname = qname
    self.seq = seq
    self.qual = qual
    self.is_read1 = is_read1
    self.is_read2 = is_read2
    self.is_paired = is_read1 or is_read2
    self.is_reverse = False
    self.is_qcfail = False
    self.is_secondary = False

  def sam_constructor(self, record):
    """Inits object using SAM record string"""
    entries = record.strip().split()
    if len(entries) < 11:
      raise Exception("ERROR: Input stream format not recognized")
    try:
      sam_flag = int(entries[1])
    except:
      raise Exception("ERROR: Input stream format not recognized")
    self.qname = entries[0]
    self.seq = entries[9]
    self.qual = entries[10]
    self.is_read1 = bool(_IS_READ1 & sam_flag)
    self.is_read2 = bool(_IS_READ2 & sam_flag)
    self.is_qcfail = bool(_IS_QC_FAIL & sam_flag)
    self.is_secondary = bool(_IS_NOT_PRIMARY_ALIGNMENT & sam_flag)
    self.is_reverse = bool(_IS_REVERSE & sam_flag)
    self.is_paired = bool(_IS_PAIRED & sam_flag)
    return self

def _reverse_complement(seq):
  """Returns reverse complement of given sequence"""
  return "".join([_COMPLEMENT[base] for base in seq[::-1]])

def _read_fastq_record(fastq_file):
  """Reads one record from a FASTQ file
  Args:
    fastq_file: Open file object
  Returns:
    record: Object with qual, seq, qname attributes.
      Will return None when EOF reached (or empty line)
  """
  qname = fastq_file.readline().strip()
  if not qname:
    return
  if qname[0] != "@":
    raise Exception("ERROR: Invalid FASTQ entry")
  else:
    qname = qname[1::]
  if qname[-1] == "1":
    is_read1 = True
    is_read2 = False
  elif qname[-1] == "2":
    is_read1 = False
    is_read2 = True
  else:
    is_read1 = False
    is_read2 = False
  seq = fastq_file.readline().strip()
  if fastq_file.readline()[0] != "+":
    raise Exception("ERROR: Invalid FASTQ entry")
  qual = fastq_file.readline().strip()
  if (qname[len(qname)-2:len(qname)] == '/1'
      or qname[len(qname)-2:len(qname)] == '/2'):
    qname = qname[0:len(qname)-2]
  return _SimpleRecord(qname, seq, qual, is_read1, is_read2)

def _open_fastq(in_path):
  """Returns compressed or uncompressed FASTQ file handle"""
  try:
    in_fq = gzopen(in_path, 'r')
    in_fq.readline()
  except IOError:
    in_fq = open(in_path, 'r')
  in_fq.seek(0)
  return in_fq

def _open_alignment_file(in_path):
  """Returns alignment file handle for BAM, SAM, or CRAM"""
  input_extension = in_path.split(".")[-1].lower()
  if input_extension == "bam":
    alignment_file = Samfile(in_path, "rb")
  elif input_extension == "sam":
    alignment_file = Samfile(in_path, "r")
  elif input_extension == "cram":
    alignment_file = Samfile(in_path, "rc")
  return alignment_file

def split_sam_stream(out_prefix, log_output=stdout, output_queue=None):
  """Reads SAM records from stdin and converts to split and compressed FASTQ
  Args:
    out_prefix: Path prefix for output files
    output_queue: Queue() object from queue builtin. If provided, output chunks
      will be added as they are created. Multithreading compatible.
  Returns:
    output_paths: List of FASTQ pair paths created
  """
  log_output.write("Splitting SAM records from stdin into paired FASTQ files "
                   "with format:\n{}_XXXXX_R[1,2].fq.gz\n".format(out_prefix))
  record_writer = _RecordWriter(out_prefix, log_output, output_queue)
  record_count = 0
  filter_secondary_count = 0
  filter_qc_count = 0
  filter_improper_pair_count = 0
  record_dict = {}
  for line in stdin:
    record = _SimpleRecord().sam_constructor(line)
    record_count += 1
    if _FILTER_NON_PRIMARY and record.is_secondary:
      filter_secondary_count += 1
    elif _FILTER_QC and record.is_qcfail:
      filter_qc_count += 1
    elif record.is_paired:
      if record.qname in record_dict:
        record2 = record_dict.pop(record.qname)
        if not ((record.is_read1 and record2.is_read2) or
                (record.is_read2 and record2.is_read1)):
          filter_improper_pair_count += 1
          log_output.write("Improper pair detected: {}\n".format(record.qname))
          log_output.write("Discarding one of these reads\n")
          record_dict[record.qname] = record
        else:
          record_writer.write_paired_records(record, record2)
      else:
        record_dict[record.qname] = record
    elif not _FILTER_UNPAIRED:
      record_writer.write_unpaired_record(record)
    if not record_count % _UPDATE_COUNT:
      log_output.write('Processed {} records\n'.format(str(record_count)))
  if not _FILTER_UNPAIRED:
    for record in record_dict.values():
      record_writer.write_unpaired_record(record)
  if _FILTER_QC:
    log_output.write("Filtered {} records failed "
                     "QC check\n".format(str(filter_qc_count)))
  if _FILTER_NON_PRIMARY:
    log_output.write("Filtered {} records that were "
                     "secondary\n".format(str(filter_secondary_count)))
  log_output.write("Filtered {} records that were improperly "
                   "paired\n".format(str(filter_improper_pair_count)))
  return record_writer.output_paths

def split_alignment_file(in_path, out_prefix, log_output=stdout,
                         output_queue=None):
  """Converts BAM, SAM, or CRAM file into split FASTQ files
  Args:
    in_path: Full path to BAM, SAM, or CRAM file to convert.
    out_prefix: Full path prefix for FASTQ files to generate.
    output_queue: Queue() object from queue builtin. If provided, output chunks
      will be added as they are created. Multithreading compatible.
  Returns:
    output_paths: List of FASTQ pair paths created
  """
  log_output.write("Input: {}\nSplitting alignment file "
                   "into paired FASTQ files with format:\n"
                   "{}_XXXXX_R[1,2].fq.gz\n".format(in_path, out_prefix))
  alignment_file = _open_alignment_file(in_path)
  record_writer = _RecordWriter(out_prefix, log_output, output_queue)
  record_count = 0
  filter_qc_count = 0
  filter_secondary_count = 0
  filter_improper_pair_count = 0
  record_dict = {}
  for record in alignment_file:
    record_count += 1
    if _FILTER_NON_PRIMARY and record.is_secondary:
      filter_secondary_count += 1
    elif _FILTER_QC and record.is_qcfail:
      filter_qc_count += 1
    elif record.is_paired:
      if record.qname in record_dict:
        record2 = record_dict.pop(record.qname)
        if not ((record.is_read1 and record2.is_read2) or
                (record.is_read2 and record2.is_read1)):
          filter_improper_pair_count += 1
          log_output.write("Improper pair detected: {}\n".format(record.qname))
          log_output.write("Discarding one of these reads\n")
          record_dict[record.qname] = record
        else:
          record_writer.write_paired_records(record, record2)
      else:
        record_dict[record.qname] = record
    elif not _FILTER_UNPAIRED:
      record_writer.write_unpaired_record(record)
    if not record_count % _UPDATE_COUNT:
      log_output.write('Processed {} records\n'.format(str(record_count)))
  if not _FILTER_UNPAIRED:
    for record in record_dict.values():
      record_writer.write_unpaired_record(record)
  if _FILTER_QC:
    log_output.write("Filtered {} records failed "
                     "QC check\n".format(str(filter_qc_count)))
  if _FILTER_NON_PRIMARY:
    log_output.write("Filtered {} records that were "
                     "secondary\n".format(str(filter_secondary_count)))
  log_output.write("Filtered {} records that were improperly "
                   "paired\n".format(str(filter_improper_pair_count)))
  return record_writer.output_paths

def split_interleaved_fastq(fq_paths, out_prefix, log_output=stdout,
                            output_queue=None):
  """Converts interleaved FASTQ files to split and paired chunks
  Args:
    fq_paths: List of paths to input interleaved FASTQ files.
    out_prefix: Prefix for output FASTQ chunks
    output_queue: Queue() object from queue builtin. If provided, output chunks
      will be added as they are created. Multithreading compatible.
  Returns:
    output_paths: List of FASTQ pair paths created
  """
  for fq_path in fq_paths:
    log_output.write("Input: {}\n".format(fq_path))
  log_output.write("Splitting interleaved FASTQ file(s) into paired FASTQ files"
                   " with format:\n{}_XXXXX_R[1,2].fq.gz\n".format(out_prefix))
  record_writer = _RecordWriter(out_prefix, log_output, output_queue)
  record_count = 0
  for fq_path in fq_paths:
    in_fq = _open_fastq(fq_path)
    while True:
      record_1 = _read_fastq_record(in_fq)
      record_2 = _read_fastq_record(in_fq)
      record_count += 1
      #EOF will return None for a record
      if not record_1 or not record_2:
        if not _FILTER_UNPAIRED:
          if record_1:
            record_writer.write_unpaired_record(record_1)
          else:
            record_writer.write_unpaired_record(record_2)
        break
      if record_1.qname != record_2.qname:
        err_msg = "Input FASTQ not sorted. Paired records not adjacent"
        log_output.write("ERROR: {}\n".format(err_msg))
        raise Exception(err_msg)
      record_writer.write_paired_records(record_1, record_2)
      if not record_count % _UPDATE_COUNT:
        log_output.write('Processed {} records\n'.format(str(record_count)))
  return record_writer.output_paths

def split_paired_fastq(fq_path_pairs, out_prefix, log_output=stdout,
                       output_queue=None):
  """Converts paired FASTQ files to split and paired chunks
  Args:
    fq_path_pairs: List of fastq pairs.
      (Ex: [["test_0_R1.fq", "test_0_R2.fq"], ["test_1_R1.fq", "test_1_R2.fq"]
    out_prefix: Prefix for output FASTQ chunks
    output_queue: Queue() object from queue builtin. If provided, output chunks
      will be added as they are created. Multithreading compatible.
  Returns:
    output_paths: List of FASTQ pair paths created
  """
  for fq_path_pair in fq_path_pairs:
    log_output.write("Input: {}\n"
                     "       {}\n".format(fq_path_pair[0], fq_path_pair[1]))
  log_output.write("Splitting paired FASTQ file(s) into smaller pairs with "
                   "format:\n{}_XXXXX_R[1,2].fq.gz\n".format(out_prefix))
  record_writer = _RecordWriter(out_prefix, log_output, output_queue)
  record_count = 0
  for fq_1_path, fq_2_path in fq_path_pairs:
    in_fq_1 = _open_fastq(fq_1_path)
    in_fq_2 = _open_fastq(fq_2_path)
    while True:
      record_1 = _read_fastq_record(in_fq_1)
      record_2 = _read_fastq_record(in_fq_2)
      record_count += 1
      #EOF will return None for a record
      if not record_1 or not record_2:
        if not _FILTER_UNPAIRED:
          while record_1:
            record_writer.write_unpaired_record(record_1)
            record_1 = _read_fastq_record(in_fq_1)
          while record_2:
            record_writer.write_unpaired_record(record_2)
            record_2 = _read_fastq_record(in_fq_2)
        break
      if record_1.qname.split()[0] != record_2.qname.split()[0]:
        err_msg = "Input FASTQ not sorted. Paired records not adjacent"
        log_output.write("ERROR: {}\n".format(err_msg))
        raise Exception(err_msg)
      record_1.is_read1 = True
      record_2.is_read2 = True
      record_writer.write_paired_records(record_1, record_2, False)
      if not record_count % _UPDATE_COUNT:
        log_output.write('Processed {} records\n'.format(str(record_count)))
  return record_writer.output_paths

def _parse_input_string(input_str, log_output):
  """Returns input string in format appropriate for splitting functions
  Function supports globbing syntax. If input is paired FASTQ files, then
  input_str should have two glob patterns separated by a comma (',').
  Because of this, paths themselves should not have commas.
  For multiple paired FASTQ files, this function assumes sorting each
  file list returned by the glob pattern will match pairs by index.
  Returns:
    String if path refers to single alignment file
    List of path strings if interleaved FASTQ file(s)
    List of path string tuples if paired FASTQ file(s)
  """
  if input_str.find(",") != -1:
    input_pattern_1, input_pattern_2 = input_str.strip().split(',')
    input_paths_1 = sorted(glob(input_pattern_1.strip()))
    input_paths_2 = sorted(glob(input_pattern_2.strip()))
    if not input_paths_1:
      err_msg = "Cannot access '{}': No such file(s)".format(input_pattern_1)
      log_output.write("ERROR: {}\n".format(err_msg))
      raise Exception(err_msg)
    if not input_paths_2:
      err_msg = "Cannot access '{}': No such file(s)".format(input_pattern_2)
      log_output.write("ERROR: {}\n".format(err_msg))
      raise Exception(err_msg)
    if len(input_paths_1) != len(input_paths_2):
      err_msg = "Globs '{}' specify different number of files".format(input_str)
      log_output.write("ERROR: {}\n".format(err_msg))
      raise Exception(err_msg)
    input_paths = []
    for index in range(len(input_paths_1)):
      input_paths.append([input_paths_1[index], input_paths_2[index]])
    return input_paths
  else:
    input_paths = glob(input_str)
    if not input_paths:
      err_msg = "Cannot access \'{}\': No such file(s)".format(input_str)
      log_output.write("ERROR: {}\n".format(err_msg))
      raise Exception(err_msg)
    if len(input_paths) == 1:
      input_path = input_paths[0]
      input_extension = input_path.split('.')[-1]
      if (input_extension == "bam" or
          input_extension == "sam" or
          input_extension == "cram"):
        return input_path
      else: #If not alignment file, then interleaved FASTQ file
        return [input_path]
    else:
      return input_paths

def split_input(input_str, out_prefix, log_output=stdout, output_queue=None):
  """Parses input string for input type (FASTQ, Alignment File) and converts
  to split and paired chunks
  Args:
    input_str: Path to input file if a single file
               Paths to input files separated by comma (',') if paired files
               Supports glob patterns for many files.
    out_prefix: Prefix for output FASTQ files
    log_output: Handle to write log output to
    output_queue: Queue() from queue builtin. If provided, finished FASTQ chunks
      will be added to the queue as they are created.
  Returns:
    output_paths: List of paired output path tuples
  """
  in_paths = _parse_input_string(input_str, log_output)
  if isinstance(in_paths, list):
    if len(in_paths[0]) == 2:
      return split_paired_fastq(in_paths, out_prefix, log_output, output_queue)
    else:
      return split_interleaved_fastq(in_paths, out_prefix,
                                     log_output, output_queue)
  else:
    return split_alignment_file(in_paths, out_prefix,
                                log_output, output_queue)

def main(input_str, out_prefix, log_path):
  log_output = open(log_path, 'w')
  start = time()
  output_paths = split_input(input_str, out_prefix, log_output)
  end = time()
  if output_paths:
    log_output.write("FASTQ splitting completed in "
                     "{} seconds\n".format(round(end-start, 2)))
    log_output.write("FASTQ pairs generated:\n{}\n".format(len(output_paths)))
  else:
    log_output.write("ERROR: No output generated\n")
  log_output.close()
  sys.exit()

if __name__ == "__main__":
  input_str = argv[1]
  out_prefix = argv[2]
  log_path = argv[3]
  main(input_str, out_prefix, log_path)
