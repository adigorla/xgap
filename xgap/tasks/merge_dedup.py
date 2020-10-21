#!/u/local/apps/python/3.7.2/bin/python3

"""Methods for merging sorted BAM files and dedupping"""


from subprocess import run, Popen, PIPE
from glob import glob
from os import path, environ, remove, fsync
from sys import stdout, argv
from sys import path as syspath
from sys import exit as sysexit
from time import time, sleep
from pysam import AlignmentFile
from tempfile import NamedTemporaryFile

#using RETRY bool to decide what INDEX_STR should be set to
try:
  XGAP_DIR = argv[9]
  RETRY = True
  syspath.insert(1,XGAP_DIR)
except IndexError:
  try:
    XGAP_DIR = argv[8]
    RETRY = False
    syspath.insert(1,XGAP_DIR)
  except IndexError:
    RETRY = False
    pass

if 'SLURM_ARRAY_TASK_ID' in  environ:
    taskid = int(environ['SLURM_ARRAY_TASK_ID'])
elif 'SGE_TASK_ID' in environ:
    if environ['SGE_TASK_ID'] != 'undefined':
        taskid = int(environ['SGE_TASK_ID'])
elif 'PBS_ARRAYID' in environ:
    taskid = int(environ['PBS_ARRAYID'])

if __name__ == "__main__":
  if RETRY:
    INDEX_STR = str(argv[7])
    JAVA_DIR = argv[8]
  else:
    INDEX_STR = str(taskid)
    JAVA_DIR = argv[7]

def split_bam(sambamaba, output_bam, out_prefix, interval_dir, log_output):
  ifiles = glob("{}/region_*.m3bp.intervals".format(interval_dir))
  exceperr = "Exception@BioD/bio/bam/reader"
  active = False
  tempf = NamedTemporaryFile(mode='w+t')

  #Using sambamba slice to split chri deduped bam
  for dfile in ifiles:
    interval = []
    f = open(dfile, 'r')
    for line in f:
      interval.append(line.strip())
    f.close
    rname = dfile.split("/")[-1].split(".")[0]
    if len(interval) > 1:
      cmdarg = [sambamba, "view", "--format=bam", "--output-filename={}_{}.bam".format(out_prefix, rname), output_bam]
      if len(interval) > 200:
        for itm in interval:
          elms = itm.strip().split(':')
          start, stop = elms.pop(-1).split('-')
          start = str(int(start)-1)
          tempf.write("{}\t{}\t{}\n".format(':'.join(elms), start, stop))
        tempf.flush()
        cmdarg.insert(-1, "-L")
        cmdarg.insert(-1, tempf.name)
      else:
        cmdarg.extend(interval)
      proc = Popen(cmdarg, stdout=PIPE, stderr=PIPE)
    else:
      proc = Popen([sambamba, "slice", "--output-filename={}_{}.bam".format(out_prefix, rname), output_bam, interval[0]], stdout=PIPE, stderr=PIPE)
    active = True
    while active == True:
      if proc.poll() is not None:
        proc.stdout.close()
        outstr = proc.stderr.read().decode('utf-8')
        if outstr.find(exceperr) == -1 and outstr.strip() != '':
          log_output.write("\n\tERROR: {}\n\n\tEcountered this unknown error while using Samabama to split BAM at {} for {}\n".format(outstr, rname, interval))
          raise Exception()
        #if no reads in region then using pysam to create a bma file with template only
        elif outstr.find(exceperr) != -1:
          infile = AlignmentFile(output_bam, "rb")
          outfile = AlignmentFile("{}_{}.bam".format(out_prefix, rname), "wb", template=infile)
          outfile.close()
          infile.close()
        proc.stderr.close()
        active = False
  tempf.close()

  return(0)

def mark_duplicates(sambamba, in_path, out_path, log_output=stdout,
                    remove_dup=False, threads=1, compress=3, hashsize=4194304, oflist=2000000, iobuffer=512):
  """Sambamba markdup
  Args:
    sambamba: Path to sambamba executable
    in_path: Path to input BAM
    out_path: Path to output deduped BAM
    log_output: File handle for log file or stdout
    remove_dup: True to remove duplicates, False otherwise.
    threads: number of threads to use
    compress: specify compression level of the resulting file (from 0 to 9)")
    hashsize: size of hash table for finding read pairs;
              will be rounded down to the nearest power of two;
              should be > (average coverage) * (insert size) for good performance
    oflist: size of overflow list where reads, thrown away from the hash table,
            get a second chance to meet their pairs. Increasing this reduces the
            number of temp files created.
    iobuffer: 2 buffer of this size (in MB) used for reading and writing BAM during
              2nd pass
  """

  tmpdir = "{}/../../xgap_intermediate_files/{}/".format(out_dir,sample_id)
  cmd = [sambamba, "markdup", "--nthreads={}".format(threads),
         "--compression-level={}".format(compress), "--tmpdir={}".format(tmpdir), "--overflow-list-size={}".format(oflist), "--io-buffer-size={}".format(iobuffer)]
  if remove_dup:
    cmd.append("--remove-duplicates")
  cmd.append(in_path)
  cmd.append(out_path)
  start = time()
  run(cmd, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("MarkDuplicates completed in {} seconds\n".format(end - start))
  log_output.flush()
  fsync(log_output.fileno())

def main(picard_jar, sample_id, n_regions, out_dir, interval_dir, log_prefix):
  # Merge bams for current region
  if ((INDEX_STR == "chri") or (int(INDEX_STR) == 0)):
    region_name = "chri"
    log_path = "{}.{}.log".format(log_prefix, region_name)
  else:
    index = int(INDEX_STR)-1
    region = str(index).zfill(len(n_regions))
    region_name = "region_{}".format(region)
    log_path = "{}.{}.log".format(log_prefix, region)
  log_output = open(log_path, 'w')
  complete_start = time()
  input_bam = "{}/Dedup/{}_{}.bam".format(out_dir, sample_id, region_name)
  input_bams = glob("{}/Bwa/{}/*.bam".format(out_dir,region_name))
  """
  start = time()
  cmd = [sambamba,"merge",
         "--nthreads=1",
         "--compression-level=3"]
  run(cmd+[input_bam]+input_bams, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("Finished merging in {} seconds\n".format(end-start))
  log_output.flush()
  fsync(log_output.fileno())
  """
  # trying to merge more than 300 file per process could trigger a Too Many Files ERROR in most systems!
  # if your system has a soft file discriptor limit >> 1024, merge all bams in one step for faster performance
  input_bams = [input_bams[x:x+300] for x in range(0, len(input_bams), 300)]
  temp_bams = ["{}/Dedup/{}_{}_{}.bam".format(out_dir, sample_id, region_name, i) for i in range(len(input_bams))]
  log_output.write("Merging BAMs for {}\n".format(region_name))
  log_output.flush()
  fsync(log_output.fileno())
  start = time()
  #using sambamaba to merge bams
  cmd = [sambamba,"merge",
         "--nthreads=1",
         "--compression-level=3"]
  #parallel merge
  processes = [Popen(cmd+[temp_bams[n]]+files , stdout=PIPE, stderr=PIPE) for n, files in enumerate(input_bams)]
  delay = 25
  while processes:
    sleep(round(delay))
    delay /= 5
    for p in processes[:]:
      if p.poll() is not None:
        log_output.write(p.stdout.read().decode('utf-8'))
        p.stdout.close()
        log_output.write(p.stderr.read().decode('utf-8'))
        p.stderr.close()
#        pend = time()
#        log_output.write("Finished {} in {} seconds\n".format(p, pend-start))
        processes.remove(p)
  step1 = time()
  log_output.write("Finished merging step1 in {} seconds\n".format(step1-start))
  log_output.flush()
  fsync(log_output.fileno())
  run(cmd+[input_bam]+temp_bams, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("Finished merging step2 in {} seconds\n".format(end-step1))
  log_output.flush()
  fsync(log_output.fileno())
  #delete temp_bams
  for f in temp_bams:
    if path.isfile(f):
      remove(f)
      remove(f+'.bai')
  log_output.write("Finished merging in {} seconds\n".format(end-start))
  log_output.flush()
  fsync(log_output.fileno())

  # Dedup merged bam
  output_bam = "{}/Dedup/{}_{}.dedup.bam".format(out_dir, sample_id,
                                                 region_name)

  log_output.write("Running MarkDuplicates on merged BAM.\n")
  log_output.flush()
  fsync(log_output.fileno())
  mark_duplicates(sambamba, input_bam, output_bam, log_output)

  # Split dedupped chri bam
  #TODO: See if you can replace current splitting code with a MUCH faster alternative
  if region_name == "chri":
    log_output.write("Splitting chri BAM\n")
    out_prefix = "{}/Dedup/{}_chri.dedup".format(out_dir, sample_id)
    split_bam(sambamba, output_bam, out_prefix, interval_dir, log_output)
  complete_end = time()
  log_output.write("Merge and dedup completed in "
                   "{} seconds\n".format(round(complete_end - complete_start)))
  log_output.flush()
  fsync(log_output.fileno())
  log_output.close()
  sysexit()

if __name__ == "__main__":
  sambamba = argv[1]
  sample_id = argv[2]
  n_regions = argv[3]
  out_dir = argv[4]
  interval_dir = argv[5]
  log_prefix = argv[6]
  main(sambamba, sample_id, n_regions, out_dir, interval_dir, log_prefix)
