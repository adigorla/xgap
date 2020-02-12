#!/usr/bin/python3

"""Apply BQSR, call variants with HaplotypeCaller and generate vcf files using GenotypeGVCF"""

from subprocess import run
from os import environ, fsync
from sys import stdout, argv
from sys import exit as sysexit
from time import time

if 'SLURM_ARRAY_TASK_ID' in  environ:
    taskid = int(environ['SLURM_ARRAY_TASK_ID'])
elif 'SGE_TASK_ID' in environ:
    taskid = int(environ['SGE_TASK_ID'])
elif 'PBS_ARRAYID' in environ:
    taskid = int(environ['PBS_ARRAYID'])

if __name__ == "__main__":
  n_regions = argv[6]
  try:
    INDEX_STR = str(argv[11])
  except IndexError:
    INDEX_STR = str(int(n_regions) - taskid)
  JAVA_DIR = argv[10]

def print_reads(gatk_jar, in_path, out_path, ref_fa, interval_path,
                n_cthreads, log_output=stdout, bqsr_table=None):
  """GATK PrintReads
  Prints reads from input BAM applying given filters or recalibration table
  Args:
    gatk_jar: Path to GATK jar file
    in_path: Input BAM path
    out_path: Output BAM path
    ref_fa: Path to reference fasta file
    interval_path: Path to interval file for input BAM
    log_output: File handle for log file.
    n_cthreads: Number of computing threads
    bqsr_table: Path to BQSR table if applying BaseRecalibrator data
  """
  cmd = [JAVA_DIR, "-Xmx5g", "-Xms512m", "-Djava.awt.headless=true", "-jar", gatk_jar,
         "-T", "PrintReads",
         "-R", ref_fa,
         "-I", in_path,
         "-L", interval_path,
         "-o", out_path,
         "-nct", n_cthreads]
  if bqsr_table:
    cmd.append("-BQSR")
    cmd.append(bqsr_table)
  start = time()
  run(cmd, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("PrintReads completed in {} seconds\n".format(end-start))
  log_output.flush()
  fsync(log_output.fileno())

def haplotype_caller(gatk_jar, in_path, out_path, ref_fa, interval_path,
                     dbsnp_path, log_output=stdout):
  """GATK HaplotypeCaller
  Calls SNPs and indels from input BAM
  Args:
    gatk_jar: Path to GATK jar file
    in_path: Path to input BAM file
    out_path: Path to output VCF file
    ref_fa: Path to reference fasta file
    interval_path: Path to interval file
    dbsnp_path: Path to dbsnp VCF file
    log_output: File handle for log file
  """
  cmd = [JAVA_DIR, "-Xmx4g", "-Xms512m", "-Djava.awt.headless=true", "-jar", gatk_jar,
         "-T", "HaplotypeCaller",
         "-R", ref_fa,
         "-I", in_path,
         "-o", out_path,
         "-L", interval_path,
         "--dbsnp", dbsnp_path,
         "-dt", "NONE",
         "-baq", "OFF",
         "-mbq", "17",
         "--excludeAnnotation", "BaseQualityRankSumTest", "--excludeAnnotation", "MappingQualityRankSumTest", "--excludeAnnotation", "ReadPosRankSumTest",
         "-stand_call_conf", "50",
         "-pairHMM", "VECTOR_LOGLESS_CACHING"]
  start = time()
  run(cmd, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("HaplotypeCaller completed in "
                   "{} seconds\n".format(end-start))
  log_output.flush()
  fsync(log_output.fileno())

def VCFgen(gatk_jar,in_gvcf,out_vcf, ref_fa,log_output):
  """GATK GenotypeGVCFs
  Performs joint genotyping on gVCF files produced by HaplotypeCaller
  Args:
    gatk_jar: Path to GATK jar file
    in_gvcF: Path to input gVCF file
    out_vcf: path to output VCF file
    ref_fa: Path to reference fasta file
    log_output: File handle for log file
  """
  cmd = [JAVA_DIR, "-Xmx4g", "-Xms512m", "-Djava.awt.headless=true", "-jar", gatk_jar,
         "-T", "GenotypeGVCFs",
         "-R", ref_fa,
         "--variant", in_gvcf,
         "--out", out_vcf,
         "--useNewAFCalculator"]
  start = time()
  run(cmd, stdout=log_output, stderr=log_output)
  end = time()
  log_output.write("GenotypeGVCFs completed in "
                   "{} seconds\n".format(end-start))
  log_output.flush()
  fsync(log_output.fileno())

def main(gatk_jar, sample_id, ref_fa, interval_dir, out_dir, n_regions,
         n_cthreads, dbsnp_path, log_prefix):
  log_path = "{}.{}.log".format(log_prefix, INDEX_STR)
  log_output = open(log_path, 'w')
  #Apply base recal
  log_output.write("Applying BQSR table\n")
  log_output.flush()
  fsync(log_output.fileno())
  bqsr_table = "{}/Recal/{}.dedup.bqsr.csv".format(out_dir, sample_id)
  region_name = INDEX_STR.zfill(len(n_regions))
  uncalibrated_bam = "{}/Recal/{}_region_{}.dedup.bam".format(out_dir, sample_id,
                                                       region_name)
  recal_bam = "{}/bam/{}_{}.dedup.recal.bam".format(out_dir, sample_id,
                                                      region_name)
  interval_path = "{}/region_{}.m2bp.intervals".format(interval_dir,
                                                       region_name)
  print_reads(gatk_jar, uncalibrated_bam, recal_bam, ref_fa, interval_path,
              n_cthreads, log_output, bqsr_table)
  #HaplotypeCaller
  log_output.write("Running HaplotypeCaller\n")
  log_output.flush()
  fsync(log_output.fileno())
  out_gvcf = "{}/vcf/{}_region_{}.g.vcf.gz".format(out_dir, sample_id, region_name)
  out_vcf = out_gvcf.replace(".g.vcf",".vcf")
  interval_path = "{}/region_{}.m1bp.intervals".format(interval_dir,
                                                       region_name)
  haplotype_caller(gatk_jar, recal_bam, out_vcf, ref_fa, interval_path,
                   dbsnp_path, log_output)
  #Genarate VCF's
  out_vcf = out_gvcf.replace(".g.vcf",".vcf")
  log_output.close()
  sysexit(0)

if __name__ == "__main__":
  gatk_jar = argv[1]
  sample_id = argv[2]
  ref_fa = argv[3]
  interval_dir = argv[4]
  out_dir = argv[5]
  n_cthreads = argv[7]
  dbsnp_path = argv[8]
  log_prefix = argv[9]
  main(gatk_jar, sample_id, ref_fa, interval_dir, out_dir, n_regions,
       n_cthreads, dbsnp_path, log_prefix)
