#!/u/local/apps/python/3.7.2/bin/python3

"""Variant Quality Score Recalibration"""


from subprocess import run
from subprocess import check_output


from os import path, environ, fsync
from sys import stdout, argv
from sys import exit as sysexit
from time import time, sleep
from pysam import Samfile, index, AlignmentFile
from glob import glob

if 'SLURM_ARRAY_TASK_ID' in  environ:
	taskid = int(environ['SLURM_ARRAY_TASK_ID'])
elif 'SGE_TASK_ID' in environ:
	if environ['SGE_TASK_ID'] != 'undefined':
		taskid = int(environ['SGE_TASK_ID'])
elif 'PBS_ARRAYID' in environ:
	taskid = int(environ['PBS_ARRAYID'])

xmx = "-Xmx10g"
xms = "-Xms8g"

if __name__ == "__main__":
	JAVA_DIR = argv[6]



def MergeVCFS(bcftools_path, in_list, out_path, log_output=stdout):
	cmd = [bcftools_path, "concat",
		"--allow-overlaps",
		"--remove-duplicates",
		"-O", "z", "-o", out_path]
	cmd.extend(in_list)
	
	start = time()
	run(cmd, stdout=log_output, stderr=log_output, check=True)
	sleep(10)
	out_idx = out_path+".tbi"
	run([bcftools_path, "index", "--force", "--tbi", "--output", out_idx, out_path], stdout=log_output, stderr=log_output)
	end = time()
	sleep(5)
	log_output.write("Merge VCFs completed in {} seconds\n".format(end-start))
	log_output.flush()
	fsync(log_output.fileno())
	sleep(1)

def IndelVariantRecalibrator(gatk_jar, in_vcf, recal_table, tranches, ref_fa, hapmap, dbsnp, omni, mills, onekindel, log_output=stdout):
	cmd = [JAVA_DIR, xms, xmx,"-Djava.awt.headless=true", "-jar", gatk_jar,
		"VariantRecalibrator",
		"-R", ref_fa,
		"-V", in_vcf,
		"-O", recal_table,
		"--tranches-file", tranches,
		"--trust-all-polymorphic",
		"-an", "DP",
		"-an", "FS",
		"-an", "MQ","-an", "QD", "-an", "SOR",
		"-mode", "INDEL",
		"-resource:hapmap,known=false,training=true,truth=true,prior=10.0", hapmap,
		"-resource:dbsnp,known=true,training=false,truth=false,prior=2.0", dbsnp,
		"-resource:omni,known=true,training=false,truth=false,prior=7.0", omni,
		"-resource:1000G,known=false,training=true,truth=false,prior=10.0", onekindel,
		"-resource:mills,known=false,training=true,truth=true,prior=12.0", mills
		]
	start = time()
	run(cmd, stdout=log_output, stderr=log_output)
	end = time()
	log_output.write("IndelVariantRecalibrator completed in {} seconds\n".format(end-start))
	log_output.flush()
	fsync(log_output.fileno())

def SNPVariantRecalibrator(gatk_jar, in_vcf, recal_table, tranches, ref_fa, hapmap, dbsnp, omni, mills, onekindel, log_output=stdout):
	cmd = [JAVA_DIR, xms, xmx,"-Djava.awt.headless=true", "-jar", gatk_jar,
		"VariantRecalibrator",
		"-R", ref_fa,
		"-V", in_vcf,
		"-O", recal_table,
		"--tranches-file", tranches,
		"--trust-all-polymorphic",
		"-an", "DP",
		"-an", "FS",
		"-an", "MQ","-an", "QD", "-an", "SOR",
		"-mode", "SNP",
		"--resource:hapmap,known=false,training=true,truth=true,prior=10.0", hapmap,
		"--resource:dbsnp,known=true,training=false,truth=false,prior=2.0", dbsnp,
		"--resource:omni,known=true,training=false,truth=false,prior=7.0", omni,
		"--resource:1000G,known=false,training=true,truth=false,prior=10.0", onekindel,
		"--resource:mills,known=false,training=true,truth=true,prior=12.0", mills
		]
	start = time()
	run(cmd, stdout=log_output, stderr=log_output)
	end = time()
	log_output.write("SNPVariantRecalibrator completed in {} seconds\n".format(end-start))
	log_output.flush()
	fsync(log_output.fileno())

def IndelApplyVQSR(gatk_jar, in_vcf, recal_table, tranches, recal_vcf, ref_fa, log_output=stdout):
	cmd = [JAVA_DIR, xms, xmx,"-Djava.awt.headless=true", "-jar", gatk_jar,
		"ApplyVQSR",
		"-R", ref_fa,
		"-V", in_vcf,
		"--recal-file", recal_table,
		"--tranches-file", tranches,
		"--truth-sensitivity-filter-level", "99.0",
		"-mode", "INDEL",
		"-O", recal_vcf]
	start = time()
	run(cmd, stdout=log_output, stderr=log_output)
	end = time()
	log_output.write("IndelApplyVQSR completed in {} seconds\n".format(end-start))
	log_output.flush()
	fsync(log_output.fileno())


def SNPApplyVQSR(gatk_jar, in_vcf, recal_table, tranches, recal_vcf, ref_fa, log_output=stdout):
	cmd = [JAVA_DIR, xms, xmx,"-Djava.awt.headless=true", "-jar", gatk_jar,
		"ApplyVQSR",
		"-R", ref_fa,
		"-V", in_vcf,
		"--recal-file", recal_table,
		"--tranches-file", tranches,
		"--truth-sensitivity-filter-level", "99.0",
		"-mode", "SNP",
		"-O", recal_vcf]
	start = time()
	run(cmd, stdout=log_output, stderr=log_output)
	end = time()
	log_output.write("SNPApplyVQSR completed in {} seconds\n".format(end-start))
	log_output.flush()
	fsync(log_output.fileno())

def main(gatk_jar, bcftools, sample_id, out_dir, log_prefix, mills, dbsnp, onekindel, hapmap, omni, ref_fa, user_option):
	log_output = open(log_prefix, 'w')
	# setting all required filenames
	vcfs = glob("{}/vcf/{}_region_[0-9]*[0-9].vcf.gz".format(out_dir, sample_id))
	merged_vcf = "{}/vcf/{}.vcf.gz".format(out_dir, sample_id)
	indel_recal_table = "{}/vcf/{}_indel.recal".format(out_dir, sample_id)
	indel_tranches = "{}/vcf/{}_indel.tranches".format(out_dir, sample_id)
	snp_recal_table = "{}/vcf/{}_snp.recal".format(out_dir, sample_id)
	snp_tranches = "{}/vcf/{}_snp.tranches".format(out_dir, sample_id)
	indel_recal_vcf = "{}/vcf/{}_recalibrated_indels.vcf.gz".format(out_dir, sample_id)
	snp_recal_vcf = "{}/vcf/{}_recalibrated_snps.vcf.gz".format(out_dir, sample_id)
	
	#Check GATK Version
	output=check_output([JAVA_DIR, "-Xmx4g", "-Xms512m","-Djava.awt.headless=true", "-jar", gatk_jar, "--version"])
	charstr=output.decode('UTF_8')
	gatk_ver="4"
	for i, c in enumerate(charstr):
		if c.isdigit():
			gatk_ver=(charstr[i])
			break
	
	#Merge VCFs
	log_output.write("Merging VCFs\n")
	log_output.flush()
	fsync(log_output.fileno())

	MergeVCFS(bcftools, vcfs, merged_vcf, log_output)
	
	#Check User Options 
	if(user_option.lower() == "n"):
		log_output.write("USER has opted out of VQSR\n")
		sysexit(0)
	elif(gatk_ver != "4"):
		log_output.write("\n\tWARNING: GATK Version incompatible with VQSR!\n\tPlease try re-running this step with GATK4.\n")
		sysexit(0) 	
	
	#Indel Variant Recal
	log_output.write("Calculating indel recal table...\n")
	log_output.flush()
	fsync(log_output.fileno())

	IndelVariantRecalibrator(gatk_jar, merged_vcf, indel_recal_table, indel_tranches, ref_fa, hapmap, dbsnp, omni, mills, onekindel,  log_output)

	#SNP Variant Recal
	log_output.write("Calculating snp recal table...\n")
	log_output.flush()
	fsync(log_output.fileno())

	SNPVariantRecalibrator(gatk_jar, merged_vcf, snp_recal_table, snp_tranches, ref_fa, hapmap, dbsnp, omni, mills, onekindel,  log_output)

	#Apply Indel VQSR
	log_output.write("Applying VQSR to Indels\n")
	log_output.flush()
	fsync(log_output.fileno())

	IndelApplyVQSR(gatk_jar, merged_vcf, indel_recal_table, indel_tranches, indel_recal_vcf, ref_fa, log_output)

	#Apply SNP VQSR
	log_output.write("Applying VQSR to SNPs\n")
	log_output.flush()
	fsync(log_output.fileno())

	SNPApplyVQSR(gatk_jar, merged_vcf, snp_recal_table, snp_tranches, snp_recal_vcf, ref_fa, log_output)

	sysexit(0)

#Need to add stop times

if __name__ == "__main__":
	gatk_jar = argv[1]
	bcftools = argv[2]
	sample_id = argv[3]
	out_dir = argv[4]
	log_prefix = argv[5]
	ref_fa = argv[7]
	onekindel = argv[8]
	hapmap = argv[9]
	omni = argv[10]
	mills = argv[11]
	dbsnp = argv[12]
	user_option = argv[13]
	main(gatk_jar, bcftools, sample_id, out_dir, log_prefix, mills, dbsnp, onekindel, hapmap, omni, ref_fa, user_option)
