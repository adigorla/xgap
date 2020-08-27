#!/ifshome/agorla/data_bucket/apps/python3.7.4/bin/python3

"""XGAP pipeline"""

# TODO: monitor
#       python binary variable

from glob import glob
from os import path, mkdir, remove
from shutil import move, rmtree
import sys
from sys import argv
from glob import glob
from sys import exit as sysexit

try:
    from xgap import scheduler
    from xgap import workflow
    try:
        config = workflow.setup_utils.read_simple_yaml(sys.path[2])
        task_scheduler = getattr(scheduler, config["scheduler"])
    except Exception:
        task_scheduler = scheduler.slurm
except ModuleNotFoundError:
    XGAP_DIR = argv[4]
    sys.path.insert(1,XGAP_DIR)
    sys.path.insert(2,argv[3])
    from xgap import scheduler
    from xgap import workflow
    config = workflow.setup_utils.read_simple_yaml(sys.path[2])
    task_scheduler = getattr(scheduler, config["scheduler"])


from xgap import tasks
from xgap.workflow.checkpoint import Task

# Directory containing task scripts
task_dir = path.dirname(tasks.__file__)

class TaskFastq(Task):
  error_terms = [
                 "ERROR",
                 "Exception"
                ]
  success_terms = [
                   "Processed",
                   "FASTQ splitting completed"
                  ]
  def run(self, rerun_indices=None):
    job_name = "fq_{}".format(self.sample_id)
    output_dir = "{}/Fastq".format(self.out_dir)
    out_prefix = "{}/{}".format(output_dir, self.sample_id)
    if not path.isdir(output_dir):
      mkdir(output_dir)
    log_dir = "{}/Fastq".format(self.log_dir)
    if not path.isdir(log_dir):
      mkdir(log_dir)
    log_output = "{}/Fastq/fq.log".format(self.log_dir)
    cmd = "{}/fastq_prep.py {} {} {}".format(task_dir, self.sample_path,
                                             out_prefix, log_output)
    mem = self.config["avail-memory"][0]
    runtime = self.config["avail-time"][0]
    num_tasks = 1
    log_output = "/dev/null"
    job_id = [task_scheduler.submit_job(job_name, cmd, mem, runtime, log_output,
                                        num_tasks)]
    submit_checkpoint(self, job_id)


  def clean_up(self):
    pass

  def gen_output_paths(self):
    # Only available after job has run (log exists)
    output_paths = [[]]
    n_chunks = self.get_num_chunks()
    for i in range(n_chunks):
      output_paths[0].append("{}/Fastq/{}_{}_{}.fq.gz".format(self.out_dir,
                                                           self.sample_id,
                                                           i, "R1"))
      output_paths[0].append("{}/Fastq/{}_{}_{}.fq.gz".format(self.out_dir,
                                                           self.sample_id,
                                                           i, "R2"))
    return output_paths

  def gen_log_paths(self):
    log_paths = [[]]
    log_paths[0].append("{}/Fastq/fq.log".format(self.log_dir))
    return log_paths

  def get_num_chunks(self):
    """Get number of FASTQ chunks from successful log file"""
    log_path = "{}/Fastq/fq.log".format(self.log_dir)
    if path.isfile(log_path):
      # Last line of log should have number of chunks
      with open(log_path, 'r') as log_file:
        for line in log_file:
          pass
      try:
          num_chunks = line.strip()
      #need to ADD ERROR HANDLING for physical hardware failure
      except NameError as e :
          reason = "\tERROR:{}\n\tYou might not have allocated enough MEMORY and/or RUN TIME for the Fastq task OR there was hardware failure. Rerun, after checking config.yml paramerters.\n".format(e) 
          if self.config['mail-notify'][0].lower() == 'y':
              msg = ["Aborted", str(self.__class__.__name__), reason, "This entire task failed!"]
              stat = super().mail_notify(msg)
              print(stat)
          print(reason)
          sysexit()
      if num_chunks.isnumeric():
        return int(num_chunks)
    return 0


class TaskBwa(Task):
  # TODO: Cleanup where qcfail and umapped moved to bam dir
  error_terms = [
                 "ERROR",
                 "Exception"
                ]
  success_terms = [
                   "Sorted",
                   "Alignment and sorting completed"
                  ]
  def run(self, rerun_indices=None):
    job_name = "bwa_{}".format(self.sample_id)
    bin_bwa = self.config["bin-bwa"]
    n_threads = self.config["n-threads"]
    ref_fa = self.config["ref-fasta"]
    working_dir = self.config["working-dir"]
    n_regions = self.config["n-regions"]
    platform_info = self.config["seq-plat"]
    interval_dir = "{}/xgap_intervals/{}/{}".format(working_dir, self.config['ref-type'], n_regions)
    log_dir = "{}/Bwa".format(self.log_dir)
    if not path.isdir(log_dir):
      mkdir(log_dir)
    log_prefix = "{}/bwa".format(log_dir)
    output_dir = "{}/Bwa".format(self.out_dir)
    if not path.isdir(output_dir):
      mkdir(output_dir)
      n_digits = len(n_regions)
      region_names = ["region_{}".format(str(i).zfill(n_digits))
                        for i in range(int(n_regions))]
      region_names.append("qcfail")
      region_names.append("unmapped")
      region_names.append("chri")
      for region_name in region_names:
        mkdir("{}/{}".format(output_dir, region_name))
    final_bam_dir = "{}/bam".format(self.out_dir)
    if not path.isdir(final_bam_dir):
      mkdir(final_bam_dir)
    mem = self.config["avail-memory"][1]
    runtime = self.config["avail-time"][1]
    log_output = "/dev/null"
    job_ids = []
    if rerun_indices:
      for index in rerun_indices:
        cmd = "{}/align_sort.py {} {} {} {} {} {} {} {} {}".format(task_dir,
                                                                bin_bwa,
                                                                n_threads,
                                                                ref_fa,
                                                                self.sample_id,
                                                                self.out_dir,
                                                                interval_dir,
                                                                log_prefix,
                                                                platform_info,
                                                                index)
        job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime,
                                                 log_output))
    else:
      cmd = "{}/align_sort.py {} {} {} {} {} {} {} {}".format(task_dir,
                                                           bin_bwa,
                                                           n_threads,
                                                           ref_fa,
                                                           self.sample_id,
                                                           self.out_dir,
                                                           interval_dir,
                                                           log_prefix,
                                                           platform_info)
      num_tasks = TaskFastq(self.sample_id, self.sample_path,
                            self.config, self.config_path).get_num_chunks()
      job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime,
                                               log_output, num_tasks))
    submit_checkpoint(self, job_ids)

  def clean_up(self):
    # Remove FASTQ dir
    fastq_dir = "{}/Fastq".format(self.out_dir)
    if path.isdir(fastq_dir):
      rmtree(fastq_dir) #If this fails, something is wrong with run()

    # Move mapq0 and unmapped to final bam dir
    bam_dir = "{}/Bwa".format(self.out_dir)
    final_bam_dir = "{}/bam".format(self.out_dir)
    for read_type in ["qcfail", "unmapped"]:
      src = "{}/{}".format(bam_dir, read_type)
      dst = "{}/{}".format(final_bam_dir, read_type)
      move(src, dst)


  def gen_output_paths(self):
    num_tasks = TaskFastq(self.sample_id, self.sample_path, self.config, self.config_path).get_num_chunks()
    n_regions = self.config["n-regions"]
    output_paths = [[] for _ in range(num_tasks)]
    n_digits = len(n_regions)
    output_dir = "{}/Bwa".format(self.out_dir)
    regions = ["qcfail", "unmapped", "chri"]
    for i in range(int(n_regions)):
      region = "region_{}".format(str(i).zfill(n_digits))
      regions.append(region)
    for region in regions:
      for index_str in range(num_tasks):
        basename = "{}_{}_R1".format(self.sample_id, index_str)
        output_path = "{}/{}/{}.bam".format(output_dir, region, basename)
        output_paths[index_str].append(output_path)
    return output_paths

  def gen_log_paths(self):
    num_tasks = TaskFastq(self.sample_id, self.sample_path, self.config, self.config_path).get_num_chunks()
    log_prefix = "{}/Bwa/bwa".format(self.log_dir)
    log_paths= []
    for index_str in range(num_tasks):
      log_paths.append(["{}.{}.log".format(log_prefix, index_str)])
    return log_paths

class TaskDedup(Task):
  error_terms = [
                 "ERROR",
                 "Exception",
                 "Unexpected",
#                 "Too many open files",
#                 "Cannot open or create file"
                ]
  success_terms = [
                   "Finished merging",
                   "MarkDuplicates completed",
                   "Merge and dedup completed"
                  ]
  def run(self, rerun_indices=None):
    job_name = "dup_{}".format(self.sample_id)
    n_regions = self.config["n-regions"]
    log_prefix = "{}/Dedup/dup".format(self.log_dir)
    log_dir = "{}/Dedup".format(self.log_dir)
    if not path.isdir(log_dir):
      mkdir(log_dir)
    output_dir = "{}/Dedup".format(self.out_dir)
    if not path.isdir(output_dir):
      mkdir(output_dir)
    sambamba = self.config["bin-sambamba"]
    working_dir = self.config["working-dir"]
    n_regions = self.config["n-regions"]
    java_dir = self.config["java-dir"]
    interval_dir = "{}/xgap_intervals/{}/{}".format(working_dir, self.config['ref-type'], n_regions)
    mem_chri = self.config["avail-memory"][2]
    runtime_chri = self.config["avail-time"][2]
    mem = self.config["avail-memory"][3]
    runtime = self.config["avail-time"][3]
    xgap_dir = self.config["xgap_path"]
    log_output = "/dev/null"
#    log_output = "/ifshome/lizhan/log/"
    job_ids = []
    if rerun_indices:
      for index in rerun_indices:
          cmd = "{}/merge_dedup.py {} {} {} {} {} {} {} {} {}".format(task_dir, 
                                                                      sambamba, 
                                                                      self.sample_id, 
                                                                      n_regions, 
                                                                      self.out_dir, 
                                                                      interval_dir, 
                                                                      log_prefix, 
                                                                      index, 
                                                                      java_dir, 
                                                                      xgap_dir)
          if index == 0 :
              job_ids.append(task_scheduler.submit_job(job_name, cmd, mem_chri, runtime_chri, log_output))
          else:
              job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime, log_output))
    else:
      cmd = "{}/merge_dedup.py {} {} {} {} {} {} {} {}".format(task_dir,
                                                         sambamba,
                                                         self.sample_id,
                                                         n_regions,
                                                         self.out_dir,
                                                         interval_dir,
                                                         log_prefix,
                                                         java_dir,
                                                         xgap_dir)
      num_tasks = int(n_regions)
      cmd_chri = "{}/merge_dedup.py {} {} {} {} {} {} {} {} {}".format(task_dir,
                                                         sambamba,
                                                         self.sample_id,
                                                         n_regions,
                                                         self.out_dir,
                                                         interval_dir,
                                                         log_prefix,
                                                         0,
                                                         java_dir,
                                                         xgap_dir)
      job_ids.append(task_scheduler.submit_job(job_name, cmd_chri, mem_chri, runtime_chri, log_output))
      job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime, log_output, num_tasks))
    submit_checkpoint(self, job_ids)
#change how args are passed
  def clean_up(self):
    # remove bwa dir
    bwa_dir = "{}/Bwa".format(self.out_dir)
    if path.isdir(bwa_dir):
      rmtree(bwa_dir)
    # remove bams that are not dedupped
    n_regions = self.config["n-regions"]
    for index in range(int(n_regions)+1):
      if index == int(n_regions):
        region_name = "chri"
      else:
        region = str(index).zfill(len(n_regions))
        region_name = "region_{}".format(region)
      output_path = "{}/Dedup/{}_{}.bam".format(self.out_dir,
                                                self.sample_id,
                                                region_name)
      try:
          remove(output_path)
      except FileNotFoundError as e:
          print("\n\t\tERROR: {} \n\t\t{} \n \t\tPipeline Execution NOT halted \n".format(e.strerror, e.filename))

  def gen_output_paths(self):
    output_paths = []
    n_regions = self.config["n-regions"]
    for index in range(int(n_regions)+1):
      output_paths.append([])
      if index == 0:
        region_name = "chri"
        for subindex in range(int(n_regions)):
          subregion_name = "region_{}".format(str(subindex).zfill(len(n_regions)))
          output_path = "{}/Dedup/{}_chri.dedup_{}.bam".format(self.out_dir,
                                                               self.sample_id,
                                                               subregion_name)
          output_paths[index].append(output_path)
      else:
        region = str(index-1).zfill(len(n_regions))
        region_name = "region_{}".format(region)
        output_path = "{}/Dedup/{}_{}.dedup.bam".format(self.out_dir,
                                                        self.sample_id,
                                                        region_name)
        output_paths[index].append(output_path)
    return output_paths

  def gen_log_paths(self):
    log_paths = []
    num_tasks = int(self.config["n-regions"]) + 1
    log_prefix = "{}/Dedup/dup".format(self.log_dir)
    for index in range(num_tasks):
        if index == 0:
            index_str = "chri"
        else:
            index_str = str(index-1).zfill(len(self.config["n-regions"]))
        log_path = "{}.{}.log".format(log_prefix, index_str)
        log_paths.append([log_path])
    return log_paths

  
class TaskRecal(Task):
  error_terms = [
                 "ERROR",
                 "Exception"
                ]
  success_terms = [ "done!",
                   "Merged chrI reads",
                   "BaseRecalibrator completed",
                   "Merge chrI and BQSR completed"
                  ]
  def run(self, rerun_indices=None):
    job_name = "rec_{}".format(self.sample_id)
    n_regions = self.config["n-regions"]
    log_prefix = "{}/Recal/rec".format(self.log_dir)
    log_dir = "{}/Recal".format(self.log_dir)
    if not path.isdir(log_dir):
      mkdir(log_dir)
    output_dir = "{}/Recal".format(self.out_dir)
    if not path.isdir(output_dir):
      mkdir(output_dir)
    gatk_jar = self.config["gatk-jar"]
    java_dir = self.config["java-dir"]
    sambamba = self.config["bin-sambamba"]
    n_threads = self.config["n-threads"]
    ref_fa = self.config["ref-fasta"]
    working_dir = self.config["working-dir"]
    interval_dir = "{}/xgap_intervals/{}/{}".format(working_dir, self.config['ref-type'], n_regions)
    known_sites = []
    known_sites.append(self.config["dbsnp-vcf"])
    known_sites.append(self.config["mills-devine-indel-vcf"])
    known_sites.append(self.config["1000g-indel-vcf"])
    known_sites_str = ",".join(known_sites)
    mem = self.config["avail-memory"][4]
    runtime = self.config["avail-time"][4]
    log_output = "/dev/null"
#    log_output = "/ifshome/lizhan/log/"
    job_ids = []
    if rerun_indices:
      for index in rerun_indices:
        cmd = "{}/base_recal.py {} {} {} {} {} {} {} {} {} {} {} {}".format(task_dir,
                                                                      gatk_jar,
                                                                      self.sample_id,
                                                                      self.out_dir,
                                                                      ref_fa,
                                                                      known_sites_str,
                                                                      interval_dir,
                                                                      n_threads,
                                                                      n_regions,
                                                                      log_prefix,
                                                                      sambamba,
                                                                      java_dir,
                                                                      index)
        job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime,
                                                 log_output))
    else:
      cmd = "{}/base_recal.py {} {} {} {} {} {} {} {} {} {} {}".format(task_dir,
                                                                 gatk_jar,
                                                                 self.sample_id,
                                                                 self.out_dir,
                                                                 ref_fa,
                                                                 known_sites_str,
                                                                 interval_dir,
                                                                 n_threads,
                                                                 n_regions,
                                                                 log_prefix,
                                                                 sambamba,
                                                                 java_dir)
      num_tasks = int(n_regions)
      job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime,
                                               log_output, num_tasks))
    submit_checkpoint(self, job_ids)

  def clean_up(self):
    #Remove Dedup dir
    dedup_dir = "{}/Dedup".format(self.out_dir)
    if path.isdir(dedup_dir):
      rmtree(dedup_dir)

  def gen_output_paths(self):
    output_paths = []
    n_regions = self.config["n-regions"]
    for index in range(int(n_regions)):
      output_paths.append([])
      region = str(index).zfill(len(n_regions))
      region_name = "region_{}".format(region)
      out_bam = "{}/Recal/{}_{}.dedup.bam".format(self.out_dir, self.sample_id,
                                                  region_name)
      out_bqsr = "{}/Recal/{}_{}.dedup.bqsr.csv".format(self.out_dir,
                                                        self.sample_id,
                                                        region)
      output_paths[index].append(out_bam)
      output_paths[index].append(out_bqsr)
    return output_paths

  def gen_log_paths(self):
    log_paths = []
    n_regions = self.config["n-regions"]
    log_prefix = "{}/Recal/rec".format(self.log_dir)
    for index in range(int(n_regions)):
      log_path = "{}.{}.log".format(log_prefix, index)
      log_paths.append([log_path])
    return log_paths

class TaskMergeRecal(Task):
  error_terms = [
                 "ERROR",
                 "Exception"
                ]
  success_terms = [
                   "GatherBqsrReports done. Elapsed time:",
                   "GatherBQSR completed"
                  ]
  def run(self, rerun_indices=None):
    job_name = "mc_{}".format(self.sample_id)
    n_regions = self.config["n-regions"]
    gatk_jar = self.config["gatk-jar"]
    java_dir = self.config["java-dir"]
    log_path = "{}/Recal/mc.log".format(self.log_dir)
    cmd = "{}/merge_recal.py {} {} {} {} {} {}".format(task_dir,
                                                    gatk_jar,
                                                    self.sample_id,
                                                    self.out_dir,
                                                    n_regions,
                                                    log_path,
                                                    java_dir)
    mem = self.config["avail-memory"][5]
    runtime = self.config["avail-time"][5]
    num_tasks = 1
    log_output = "/dev/null"
    job_id = [task_scheduler.submit_job(job_name, cmd, mem, runtime, log_output,
                                        num_tasks)]
    submit_checkpoint(self, job_id)


  def clean_up(self):
    # Remove regional BQSR Tables
    n_regions = self.config["n-regions"]
    for index in range(int(n_regions)):
      region_name = str(index).zfill(len(n_regions))
      bqsr_path = "{}/Recal/{}_{}.dedup.bqsr.csv".format(self.out_dir, self.sample_id,
                                                         region_name)
      if path.isfile(bqsr_path):
        remove(bqsr_path)

  def gen_output_paths(self):
    return [["{}/Recal/{}.dedup.bqsr.csv".format(self.out_dir, self.sample_id)]]

  def gen_log_paths(self):
    return [["{}/Recal/mc.log".format(self.log_dir)]]

class TaskBqsrHc(Task):
  error_terms = [
                 "ERROR",
                 "Exception",
                 "USER ERROR"
                ]
  success_terms = [
                   "done",
                   "PrintReads completed in",
                   "Done.",
                   "HaplotypeCaller completed in",
                   "GenotypeGVCFs completed"
                  ]
  def run(self, rerun_indices=None):
    job_name = "hc_{}".format(self.sample_id)
    n_regions = self.config["n-regions"]
    gatk_jar = self.config["gatk-jar"]
    java_dir = self.config["java-dir"]
    log_prefix = "{}/BqsrHc/hc".format(self.log_dir)
    log_dir = "{}/BqsrHc".format(self.log_dir)
    if not path.isdir(log_dir):
      mkdir(log_dir)
    output_dir = "{}/vcf".format(self.out_dir)
    if not path.isdir(output_dir):
      mkdir(output_dir)
    n_threads = self.config["n-threads"]
    ref_fa = self.config["ref-fasta"]
    working_dir = self.config["working-dir"]
    interval_dir = "{}/xgap_intervals/{}/{}".format(working_dir, self.config['ref-type'], n_regions)
    dbsnp_path = self.config["dbsnp-vcf"]
    mem = self.config["avail-memory"][6]
    runtime = self.config["avail-time"][6]
    log_output = "/dev/null"
    job_ids = []
    if rerun_indices:
      for index in rerun_indices:
        cmd = "{}/bqsr_hc.py {} {} {} {} {} {} {} {} {} {} {}".format(task_dir,
                                                                   gatk_jar,
                                                                   self.sample_id,
                                                                   ref_fa,
                                                                   interval_dir,
                                                                   self.out_dir,
                                                                   n_regions,
                                                                   n_threads,
                                                                   dbsnp_path,
                                                                   log_prefix,
                                                                   java_dir,
                                                                   index)
        job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime,
                                                 log_output))
    else:
      cmd = "{}/bqsr_hc_vcf.py {} {} {} {} {} {} {} {} {} {}".format(task_dir,
                                                              gatk_jar,
                                                              self.sample_id,
                                                              ref_fa,
                                                              interval_dir,
                                                              self.out_dir,
                                                              n_regions,
                                                              n_threads,
                                                              dbsnp_path,
                                                              log_prefix,
                                                              java_dir)
      num_tasks = int(n_regions)
      job_ids.append(task_scheduler.submit_job(job_name, cmd, mem, runtime,
                                               log_output, num_tasks))
    submit_checkpoint(self, job_ids)

  def clean_up(self):
      recal_dir = "{}/Recal".format(self.out_dir)
      if path.isdir(recal_dir):
          rmtree(recal_dir)
  
  def gen_output_paths(self):
    output_paths = []
    n_regions = self.config["n-regions"]
    for index in range(int(n_regions)):
      output_paths.append([])
      region_name = str(index).zfill(len(n_regions))
      out_bam = "{}/bam/{}_{}.dedup.recal.bam".format(self.out_dir, self.sample_id,
                                                      region_name)
      out_vcf = "{}/vcf/{}_region_{}.vcf.gz".format(self.out_dir,
                                                    self.sample_id,
                                                    region_name)
      output_paths[index].append(out_bam)
      output_paths[index].append(out_vcf)
    return output_paths

  def gen_log_paths(self):
    log_paths = []
    n_regions = self.config["n-regions"]
    log_prefix = "{}/BqsrHc/hc".format(self.log_dir)
    for index in range(int(n_regions)):
      log_path = "{}.{}.log".format(log_prefix, index)
      log_paths.append([log_path])
    return log_paths

class TaskVQSR(Task):
  error_terms = { 
                  "ERROR",
                  "Exception"
                }
  success_terms = {}
  def run(self, rerun_indicies=None):
    job_name="VQSR_{}".format(self.sample_id)
    log_path = "{}/BqsrHC/vqsr.log".format(self.log_dir)
    log_dir="{}/BqsrHC/".format(self.log_dir)
    if not path.isdir(log_dir):
      mkdir(log_dir)
    working_dir = self.config["working-dir"]
    dbsnp = self.config["dbsnp-vcf"]
    1000g = self.config["1000g-indel-vcf"]
    mills = self.config["mills-devine-indel-vcf"]
    gatk_jar=self.config["gatk-jar"]
    java_dir = self.config["java-dir"]
    bcftools_path=self.config["bcftools-path"]
    cmd = "{}/VQSR.py {} {} {} {} {}".format(task_dir,
                                                   gatk_jar,
                                                   bcftools_path,
                                                   self.sample_id,
                                                   self.out_dir,
                                                   log_path, mills, dbsnp, 1000g)
    mem = self.config["avail-memory"][7]
    runtime = self.config["avail-time"][7]
    num_tasks = 1
    log_output = "/dev/null"
    job_ids = [task_scheduler.submit_job(job_name, cmd, mem, runtime, log_output,
                                         num_tasks)]
    submit_checkpoint(self, job_ids)
  def clean_up(self):
    pass
 
  def gen_output_paths(self):
    output_paths = [[]]
    output_paths[0].append("{}/vcf/{}.vcf.gz".format(self.out_dir, self.sample_id))
    output_paths[0].append("{}/vcf/{}_excesshet.vcf.gz".format(self.out_dir, self.sample_id))
    output_paths[0].append("{}/vcf/{}_sitesonly.vcf.gz".format(self.out_dir, self.sample_id))
    output_paths[0].append("{}/vcf/{}_indels.recal".format(self.out_dir, self.sample_id))
    output_paths[0].append("{}/vcf/{}_indels.tranches".format(self.out_dir, self.sample_id))
    output_paths[0].append("{}/vcf/{}_snps.tranches".format(self.out_dir, self.sample_id))
    output_paths[0].append("{}/vcf/{}_snps.recal".format(self.out_dir, self.sample_id))
    #Append all output paths, all appended to output_paths[0]
    output_paths[0].append("{}/vcf/{}_snp.recalibrated.vcf.gz".format(self.out_dir,
                                                           self.sample_id))
    output_paths[0].append("{}/vcf/{}_indel.recalibrated.vcf.gz".format(self.out_dir,
                                                           self.sample_id))
    return output_paths

  def gen_log_paths(self):
    log_paths = [["{}/BqsrHC/vqsr.log".format(self.log_dir)]]
    return log_paths	


class TaskMergeBams(Task):
  error_terms = [
                 "ERROR",
                 "Exception"
                ]
  success_terms = [
                   "Merged regional",
                   "Merged mapq0",
                   "Merged unmapped",
                   "Merged all"
                  ]
  def run(self, rerun_indices=None):
    job_name = "mb_{}".format(self.sample_id)
    n_regions = self.config["n-regions"]
    log_path = "{}/MergeBams/mb.log".format(self.log_dir)
    log_dir = "{}/MergeBams".format(self.log_dir)
    if not path.isdir(log_dir):
      mkdir(log_dir)
    working_dir = self.config["working-dir"]
    interval_dir = "{}/xgap_intervals/{}/{}".format(working_dir, self.config['ref-type'], n_regions)
    cmd = "{}/merge_bams.py {} {} {} {} {}".format(task_dir,
                                                   self.sample_id,
                                                   self.out_dir,
                                                   n_regions,
                                                   interval_dir,
                                                   log_path)
    mem = self.config["avail-memory"][7]
    runtime = self.config["avail-time"][7]
    num_tasks = 1
    log_output = "/dev/null"
    job_ids = [task_scheduler.submit_job(job_name, cmd, mem, runtime, log_output,
                                         num_tasks)]
    submit_checkpoint(self, job_ids)


  def clean_up(self):
      
    #Successful so remove logs, commented out for testing

    #for step in ["Recal", "BqsrHc", "MergeBams", "Dedup", "Bwa", "Fastq"]:
    #  rmtree("{}/{}".format(self.log_dir, step))
     
      # remove all regional, mapq0 and unmapped bams and index files
    
    for index in range(int(self.config["n-regions"])):
        bam_path = "{}/bam/{}_{}.dedup.recal.bam".format(self.out_dir, self.sample_id, str(index).zfill(len(self.config["n-regions"])))
        bai_path = "{}/bam/{}_{}.dedup.recal.bai".format(self.out_dir, self.sample_id, str(index).zfill(len(self.config["n-regions"])))
        if path.isfile(bam_path):
            remove(bam_path)
        if path.isfile(bai_path):
            remove(bai_path)
    if path.isdir("{}/bam/qcfail".format(self.out_dir)):
        rmtree("{}/bam/qcfail".format(self.out_dir))
    if path.isdir("{}/bam/unmapped".format(self.out_dir)):
         rmtree("{}/bam/unmapped".format(self.out_dir))
             
    

  def gen_output_paths(self):
    output_paths = [[]]
    output_paths[0].append("{}/bam/{}.dedup.recal.bam".format(self.out_dir,
                                                           self.sample_id))
    output_paths[0].append("{}/bam/{}_qcfail.bam".format(self.out_dir,
                                                       self.sample_id))
    output_paths[0].append("{}/bam/{}_unmapped.bam".format(self.out_dir,
                                                                 self.sample_id))
    return output_paths

  def gen_log_paths(self):
    log_paths = [["{}/MergeBams/mb.log".format(self.log_dir)]]
    return log_paths

class TaskFinish(Task):
  error_terms = []
  success_terms = []

  def run(self):
      self.clean_up()
      flag_list = glob("{}/task_*_begin.flag".format(self.log_dir))
      flag_list = sorted(flag_list)
      remove(flag_list[-1])
      finish_flag = "{}/pipeline_finish.flag".format(self.log_dir)
      msg = ["Completed Successfully", None, None, None]
      statement = "{} Analysis {}!\n".format(self.sample_id, msg[0])
      with open(finish_flag, 'w') as finish_file:
          finish_file.write(statement)
          stat = super().mail_notify(msg)
          finish_file.write("\n{}\n".format(stat))
      print(statement)
      print(self.__class__.__name__)
      sysexit()

  def clean_up(self):
      
      #remove all bams during testing
      #CHANGE THIS 
      bam_dir = "{}/bam".format(self.out_dir)
      if path.isdir(bam_dir):
          rmtree(bam_dir)
      
    #Successful so remove logs, commented out for testing

    #for step in ["Recal", "BqsrHc", "MergeBams", "Dedup", "Bwa", "Fastq"]:
    #  rmtree("{}/{}".format(self.log_dir, step))  

  def gen_output_paths(self):
    output_paths = [[]]
    return output_paths

  def gen_log_paths(self):
    log_paths = [[]]
    return log_paths

def submit_checkpoint(task, hold_ids):
  """Submits this script as a job that runs after given hold_id"""
  job_name = "check_{}".format(task.sample_id)
  script = path.abspath(__file__)
  cmd = " ".join([script, task.sample_id, task.sample_path, task.config_path, task.xgap_path])
  mem = task.config["avail-memory"][8]
  runtime = task.config["avail-time"][8]
  log_output = task.log_dir
  # TODO: Abstract for other schedulers
  hold_ids_str = ":".join(hold_id.strip() for hold_id in hold_ids)
  task_scheduler.submit_job(job_name, cmd, mem, runtime, log_output,
                            hold_id=hold_ids_str)

def main(sampleid, samplepath,configpath):
    pipeline = [
                TaskFastq,
                TaskBwa,
                TaskDedup,
                TaskRecal,
                TaskMergeRecal,
                TaskBqsrHc,
		TaskVQSR,
#                TaskMergeBams,
                TaskFinish
               ]
    sample_id = sampleid
    sample_path = samplepath
    config_path = sys.path[2]
    config = workflow.setup_utils.read_simple_yaml(config_path)
    config["xgap_path"] = sys.path[1]
    global task_scheduler 
    task_scheduler = getattr(scheduler, config["scheduler"])
    workflow.checkpoint.main(sample_id, sample_path, config, config_path, pipeline)

if __name__ == "__main__":
  pipeline = [
              TaskFastq,
              TaskBwa,
              TaskDedup,
              TaskRecal,
              TaskMergeRecal,
              TaskBqsrHc,
	      TaskVQSR,
#              TaskMergeBams,
              TaskFinish
              ]
  sample_id = argv[1]
  sample_path = argv[2]
  config_path = sys.path[2]
  config["xgap_path"] = sys.path[1]
  workflow.checkpoint.main(sample_id, sample_path, config, config_path, pipeline)
