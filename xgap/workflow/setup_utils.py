"""Methods for setting up xgap and reading configuration settings"""

import os

__all__ = ['read_batch_file', 'read_simple_yaml', 'overwrite_config',
           'get_sequences', 'generate_intervals', 'setup_directories']

def setup_directories(config, samples, mail):
  """Creates all necessary directories for pipeline

  Possible change: Take actual directory strings, instead of depending on
    config keys.
  TODO: Improve security of mail config data

  Args:
    config: dict produced by read_simple_yaml
    samples: list of sample tuples - (sample_id, sample_path)
    mail: The configrations for email notifications
  """
  # make output dir
  output_dir = "{}/xgap_output/".format(config["output-dir"])
  if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
  # make working dir
  working_dir = "{}/xgap_intermediate_files/".format(config["working-dir"])
  if not os.path.isdir(working_dir):
    os.makedirs(working_dir)
  # make log dir
  log_dir = "{}/xgap_logs/".format(config["log-dir"])
  if not os.path.isdir(log_dir):
    os.makedirs(log_dir)
  for sample_id, sample_path in samples:
    sample_log_dir = "/".join([log_dir, sample_id])
    sample_working_dir = "/".join([working_dir, sample_id])
    sample_output_dir = "/".join([output_dir, sample_id])
    if not os.path.isdir(sample_log_dir):
      os.mkdir(sample_log_dir)
    if not os.path.isdir(sample_working_dir):
      os.mkdir(sample_working_dir)
    if not os.path.isdir(sample_output_dir):
      os.mkdir(sample_output_dir)
    if mail[-1] != None:
      with open("{}/mail.dat".format(sample_log_dir), 'w') as mfile:
        for line in mail:
          mfile.write("{}\n".format(line))
      os.chmod("{}/mail.dat".format(sample_log_dir), 0o600)
      

def read_batch_file(batch_path):
  """Reads batch file, returning list of sample ID and sample path tuples

  Args:
    batch_path: full path to batch file
  Returns:
    samples: List of tuples ('sample_id', 'sample_path')
  """
  samples = []
  with open(batch_path, 'r') as batch_file:
    for line_num, line in enumerate(batch_file):
      try:
        sample_id, sample_path = line.strip().split()
      except:
        raise Exception("Error in line {}. Expecting 2 values: "
                        "Sample_ID Sample_path".format(line_num))
      samples.append([sample_id, sample_path])
    return samples

def read_simple_yaml(yaml_path):
  """Reads yaml configuration file. Only supports simple key value assignments

  Reading key-value assignments in yaml like format without external
  dependencies.

  Args:
    yaml_path: Path to configuration file
  Returns:
    config: Dict with key value pairs from yaml file
  """
  with open(yaml_path, 'r') as yaml_file:
    config = {}
    for line in yaml_file:
      line = line.strip()
      if len(line) < 3 or line[0] == '#' or line == "---":
        continue
      if line == "...":
        break
      key, value = line.split(":", maxsplit=1)
      key = key.strip()
      config[key] = value.strip().strip("\"")
      if config[key][0] == "[":
        list_val = (("".join(config[key].strip("[").strip("]").replace('"', '').split())).split(","))
        config[key] = list_val
    return config

def overwrite_config(yaml_path, config):
  """Overwrites yaml configuration file using simple config dict

  Iterates through lines of existing config yml file. If key of current line
  is in 'config', then update value to config[key]. Otherwise, keep existing
  value.

  Args:
    yaml_path: Path to output configuration file
    config: Dictionary with key value pairs to write
  """
  temp_yaml_path = '{}.temp'.format(yaml_path)
  with open(yaml_path, 'r') as in_file, open(temp_yaml_path, 'w') as out_file:
    for line in in_file:
      line = line.strip()
      if len(line) < 3 or line[0] == '#' or line == "---" or line == "...":
        out_file.write('{}\n'.format(line))
      else:
        key = line.split(":")[0].strip()
        if key in config:
          out_file.write('{}: \"{}\"\n'.format(key, config.pop(key)))
        else:
          out_file.write('{}\n'.format(line))
  os.rename(temp_yaml_path, yaml_path)

def get_sequences(reference_dict_path):
  """Retrieves chromosome names and lengths from reference dictionary

  Args:
    reference_dict_path: Path to reference.dict
  Returns:
    sequence_ids: A list of sequence names in order given in reference dict
    sequence_lens: A python dict with sequence ids as keys and their respective
      lengths as values.
  """
  sequence_ids = []
  sequence_lens = {}
  with open(reference_dict_path, 'r') as dict_file:
    for line in dict_file:
      entries = line.strip().split()
      if entries[0] == "@SQ":
        for entry in entries[1::]:
          secondary_entries = entry.split(':', maxsplit=1)
          if secondary_entries[0].strip() == "SN":
            sequence_name = secondary_entries[1].strip()
          if secondary_entries[0].strip() == "LN":
            sequence_len = int(secondary_entries[1].strip())
        try:
          sequence_ids.append(sequence_name)
          sequence_lens[sequence_name] = sequence_len
        except:
          raise Exception("Ref dict does not have properly formatted SQ line")
  return sequence_ids, sequence_lens

def generate_intervals(reference_dict_path, interval_dir, n_regions):
  """Define equally spaced intervals from reference genome and write to files

  Will generate 4 files for each region:
    region_XXXX.intervals = Non-overlapping regions
    region_XXXX.m1bp.intervals = Regions with 1kb overlap
    region_XXXX.m2bp.intervals = Regions with 2kb overlap
    region_XXXX.m3bp.intervals = Regions with 3kb overlap
  Where 'XXXX' is the region index. The number of digits will be equal to the
  length of the string representation of n_regions.

  Args:
    reference_dict: Full path to ref.dict file containing sequence lengths
    interval_dir: Directory to write interval files.
    n_regions: Number of regions to generate
  """
  n_regions = int(n_regions)
  if n_regions <= 0:
    raise Exception("Specified number of regions not a positive integer")
  num_digits = len(str(n_regions))
  sequence_ids, sequence_lens = get_sequences(reference_dict_path)
  reference_length = sum(sequence_lens.values())
  interval_length = -(-reference_length // n_regions)
  curr_chr = sequence_ids[0]
  curr_pos = 0
  for interval in range(n_regions):
    sub_intervals = []
    finished = False
    remaining_interval = interval_length
    next_pos = curr_pos + interval_length
    while next_pos > sequence_lens[curr_chr]:
      sub_intervals.append([curr_chr, (curr_pos + 1), sequence_lens[curr_chr]])
      remaining_interval -= (sequence_lens[curr_chr] - curr_pos)
      next_chr_index = sequence_ids.index(curr_chr) + 1
      curr_pos = 0
      if remaining_interval <= 0 or next_chr_index == len(sequence_ids):
        finished = True
        break
      curr_chr = sequence_ids[next_chr_index]
      next_pos = remaining_interval
    if not finished:
      sub_intervals.append([curr_chr, (curr_pos + 1), next_pos])
      curr_pos = next_pos
    index = str(interval).zfill(num_digits)
    intv_path = '{}/region_{}.m0bp.intervals'.format(interval_dir, index)
    m1bp_path = '{}/region_{}.m1bp.intervals'.format(interval_dir, index)
    m2bp_path = '{}/region_{}.m2bp.intervals'.format(interval_dir, index)
    m3bp_path = '{}/region_{}.m3bp.intervals'.format(interval_dir, index)
    with open(intv_path, 'w') as intv_file, \
         open(m1bp_path, 'w') as m1bp_file, \
         open(m2bp_path, 'w') as m2bp_file, \
         open(m3bp_path, 'w') as m3bp_file:
      for sub_interval in sub_intervals:
        intv_file.write("{s_i[0]}:{s_i[1]}-{s_i[2]}\n".format(s_i=sub_interval))
      for sub_interval in _overlap_intervals(sub_intervals, 1000,
                                             sequence_ids, sequence_lens):
        m1bp_file.write("{s_i[0]}:{s_i[1]}-{s_i[2]}\n".format(s_i=sub_interval))
      for sub_interval in _overlap_intervals(sub_intervals, 2000,
                                             sequence_ids, sequence_lens):
        m2bp_file.write("{s_i[0]}:{s_i[1]}-{s_i[2]}\n".format(s_i=sub_interval))
      for sub_interval in _overlap_intervals(sub_intervals, 3000,
                                             sequence_ids, sequence_lens):
        m3bp_file.write("{s_i[0]}:{s_i[1]}-{s_i[2]}\n".format(s_i=sub_interval))

def _overlap_intervals(sub_intervals, overlap, sequence_ids, sequence_lens):
  """Helper function for 'generate_intervals'. Generates interval with range
  increased by given overlap on both ends

  Args:
    sub_intervals: List of lists. Each sublist defines a subinterval in the
      following form: [chr_id, start, end]
    overlap: Number of bases to increase range of sub_intervals by on each
      side. For example: ['chr1', 10, 12] => ['chr1', 8, 14] if overlap=2
    sequence_ids: List of sequence names in order
    sequence_lens: Dict of sequence lengths as values for sequence_id keys
  Returns:
    overlap_sub_intervals: list with same structure as sub_intervals with
      total range increased by overlap on each end.
  """
  overlap_sub_intervals = []
  for curr_index, sub_interval in enumerate(sub_intervals):
    if curr_index != 0 and curr_index != (len(sub_intervals) - 1):
      overlap_sub_intervals.append(sub_interval)
    else:
      new_start_chr = sub_interval[0]
      new_start_pos = sub_interval[1]
      if curr_index == 0:
        remaining_overlap = overlap
        new_start_pos -= remaining_overlap
        while new_start_pos < 1:
          prev_chr_index = sequence_ids.index(new_start_chr) - 1
          if prev_chr_index == -1:
            new_start_pos = 1
            break
          remaining_overlap = abs(new_start_pos) + 1
          new_start_chr = sequence_ids[prev_chr_index]
          new_start_pos = sequence_lens[new_start_chr] - remaining_overlap + 1
      new_end_chr = sub_interval[0]
      new_end_pos = sub_interval[2]
      if curr_index == (len(sub_intervals) - 1):
        remaining_overlap = overlap
        new_end_pos += overlap
        while new_end_pos > sequence_lens[new_end_chr]:
          next_chr_index = sequence_ids.index(new_end_chr) + 1
          if next_chr_index == len(sequence_ids):
            new_end_pos = sequence_lens[new_end_chr]
            break
          remaining_overlap = new_end_pos - sequence_lens[new_end_chr]
          new_end_chr = sequence_ids[next_chr_index]
          new_end_pos = remaining_overlap
      if new_start_chr != new_end_chr:
        overlap_sub_intervals.append([new_start_chr, new_start_pos,
                                      sequence_lens[new_start_chr]])
        for seq_index in range(sequence_ids.index(new_start_chr) + 1,
                               sequence_ids.index(new_end_chr)):
          overlap_sub_intervals.append([sequence_ids[seq_index], 1,
                                        sequence_lens[sequence_ids[seq_index]]])
        overlap_sub_intervals.append([new_end_chr, 1, new_end_pos])
      else:
        overlap_sub_intervals.append([new_end_chr, new_start_pos, new_end_pos])
  return overlap_sub_intervals
