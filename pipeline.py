#!/usr/bin/python

import sys, getopt, os
from subprocess import call

class Sample(object):
   name = ''
   tumor = ''
   normal = ''
   directory = ''

   def __init__(self, name, tumor, normal, directory):
      self.name = name
      self.tumor = tumor
      self.normal = normal
      self.directory = directory

def create_sample(line, directory):
      l = line.split()
      if l[3] == 'exome' and l[4] == 'fastq':
         sample = Sample(l[0], l[1], l[2], directory)
         return sample
      else:
         raise NotSampleError() #if the line isn't a sample we are going to use either because it is wgs or already in bam format or in the wrong format, don't make an object

def write_pbs_file(output_dir, sample, genome, username, current_dir): #TODO add timing commands to this
   file = open(output_dir + '/' + sample.name + '_pipeline.pbs')

   #default .pbs job header
   file.write('#!/bin/bash')
   file.write('#PBS -q hotel')
   file.write('#PBS -l nodes=1:ppn=8')
   file.write('#PBS -l walltime=20:00:00')
   file.write('#PBS -m bea')
   file.write('#PBS -M ' + username + '@ucsd.edu')  #assuming for now that email is of format username@ucsd.edu
   file.write('#PBS -o ' + sample.name + '_pipeline.o')
   file.write('#PBS -e ' + sample.name + '_pipeline.e')
   file.write('#PBS -N ' + sample.name + '_pipeline_' + username)
   file.write('#PBS -V')
   file.write('')

   #need to load R in the job for DNAcopy to work
   file.write('module load R')
   file.write('')

   #fastqc
   file.write('fastqc' + sample.directory + '/' + sample.tumor + '_1.fastq.gz -o ' + output_dir)
   file.write('fastqc' + sample.directory + '/' + sample.tumor + '_2.fastq.gz -o ' + output_dir)
   file.write('fastqc' + sample.directory + '/' + sample.normal + '_1.fastq.gz -o ' + output_dir)
   file.write('fastqc' + sample.directory + '/' + sample.normal + '_2.fastq.gz -o ' + output_dir)
   file.write('')

   #TODO add QC check

   #bwa for tumor
   file.write('bwa mem -t 8 -k 35 ' + genome + ' ' + sample.directory + '/' + sample.tumor + '_1.fastq.gz ' + sample.directory + '/' + sample.tumor + '_2.fastq.gz | samtools sort -@ 8 -o ' + output_dir + '/' + sample.name + '_tumor.bam -')
   file.write('/opt/biotools/GenomeAnalysisTK/4.0.11.0/gatk MarkDuplicates -I ' + output_dir + '/' + sample.name + '_tumor.bam -O ' + output_dir + '/' + sample.name + '_tumor_md.bam -M ' + output_dir + '/' + sample.name + '_tumor_md.metrics --REMOVE_DUPLICATES')
   file.write('samtools index ' + output_dir + '/' + sample.name + '_tumor_md.bam')
   file.write('')

   #bwa for normal
   file.write('bwa mem -t 8 -k 35 ' + genome + ' ' + sample.directory + '/' + sample.normal + '_1.fastq.gz ' + sample.directory + '/' + sample.normal + '_2.fastq.gz | samtools sort -@ 8 -o ' + output_dir + '/' + sample.name + '_normal.bam -')
   file.write('/opt/biotools/GenomeAnalysisTK/4.0.11.0/gatk MarkDuplicates -I ' + output_dir + '/' + sample.name + '_normal.bam -O ' + output_dir + '/' + sample.name + '_normal_md.bam -M ' + output_dir + '/' + sample.name + '_normal_md.metrics --REMOVE_DUPLICATES')
   file.write('samtools index ' + output_dir + '/' + sample.name + '_normal_md.bam')
   file.write('')

   #mpileups
   file.write('samtools mpileup -f ' + genome + ' -q 5 -Q 0 -d 50000 -B -o ' + output_dir + '/' + sample.name + '_tumor.vcf ' + output_dir + '/' + sample.name + '_tumor_md.bam')
   file.write('samtools mpileup -f ' + genome + ' -q 5 -Q 0 -d 50000 -B -o ' + output_dir + '/' + sample.name + '_normal.vcf ' + output_dir + '/' + sample.name + '_normal_md.bam')
   file.write('')

   #generate large segments for amplification and deletion calling
   file.write('varscan copynumber ' + output_dir + '/' + sample.name + '_normal.vcf ' + output_dir + '/' + sample.name + '_tumor.vcf ' + output_dir + '/' + sample.name + '_varscan --min-base-qual 0 --min-map-qual 0 --min-segment-size 50 --max-segment-size 1000')
   file.write('varscan copyCaller ' + output_dir + '/' + sample.name + '_varscan.copynumber --output-file ' + output_dir + '/' + sample.name + '_varscan.copynumber.called')
   file.write('Rscript --vanilla ' + current_dir + '/DNAcopy.R ' + output_dir + '/' + sample.name + '_varscan.copynumber.called')
   file.write('')
   
   #generate short segments for plotting so we don't have to plot every single point
   file.write('varscan copynumber ' + output_dir + '/' + sample.name + '_normal.vcf ' + output_dir + '/' + sample.name + '_tumor.vcf ' + output_dir + '/' + sample.name + '_varscan_plotting --min-base-qual 0 --min-map-qual 0 --min-segment-size 5 --max-segment-size 50')
   file.write('varscan copyCaller ' + output_dir + '/' + sample.name + '_varscan_plotting.copynumber --output-file ' + output_dir + '/' + sample.name + '_varscan_plotting.copynumber.called')
   file.write('Rscript --vanilla ' + current_dir + '/DNAcopy.R ' + output_dir + '/' + sample.name + '_varscan_plotting.copynumber.called')
   file.write('')

   #generate a plot of the sample highlighting regions of amplification and deletion
   #TODO write this python file to make a plot using numpy similar to how we've made plots in matlab
   file.write('python sample_visualization.py ' + output_dir + '/' + sample.name + '_varscan.copynumber.called.segmentation ' + output_dir + '/' + sample.name + '_varscan_plotting.copynumber.called.segmentation')
   file.write('')

   file.close()


def main(argv):
   username = os.environ.get('USER')
   current_dir = os.getcwd()

   #process inputs
   inputfile = ''
   outputfile = ''
   genome = ''
   try:
      opts, args = getopt.getopt(argv,'hi:o:g:',['help','input=','output=','genome='])
   except getopt.GetoptError:
      print('pipeline.py -i <inputfile> -o <outputfile> -g <referencegenome>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h', '--help'):
         print('pipeline.py -i <inputfile> -o <outputfile> -g <referencegenome>')
         print('pipeline.py runs the copy number pipeline developed by the BENG187 team under Dr. Alexandrov in 2019-20.')
         print('options are:')
         print('  -h or --help')
         print('  -i or --input')
         print('     The filename of the table for paired tumor and normal samples in Phoebe\'s format.')
         print('  -o or --output')
         print('     The folder name where each paired tumor-normal sample analysis files will be stored.')
         print('     Default is the name of the folder containing the input file.')
         print('  -g or --genome')
         print('     The .fa reference genome file.')
         sys.exit()
      elif opt in ('-i', '--input'):
         inputfile = arg
      elif opt in ('-o', '--output'):
         outputfile = arg
      elif opt in ('-g', '--genome'):
         genome = arg
   if inputfile == '':
      print('No input file provided')
      sys.exit(2)
   inputfolder = inputfile.rsplit('/', 1)[0] #the folder with the inputfile
   if outputfile == '':
      outputfile = '/oasis/tscc/scratch/' + username + '/' + inputfolder.rsplit('/', 1)[1] #the name of the samples folder is used as the defaul for the name of the output folder
      print('No output file provided, will output to ' + outputfile)
   if genome == ''
      print('No reference genome provided')
      sys.exit(2)
   elif genome.rsplit('.', 1)[1] != 'fa':
      print('Reference genome is not a .fa file')
      sys.exit(2)

   #make output folder
   call(['mkdir', '-p', outputfile]) #-p allows mkdir to not return an error if there is already a file named this

   #iterate over input file to get list of samples
   sample_list = []
   with open(inputfile) as f:
      for line in f:
         try:
            sample_list += create_sample(line.rstrip('\n'), inputfolder + '/fastq_files')
         except NotSampleError:
            continue

   for sample in sample_list:
      sample_dir = outputfile + '/' + sample.name
      call(['mkdir', '-p', sample_dir])
      call(['mkdir', '-p', sample_dir + '/source'])
      write_pbs_file(sample_dir, sample, genome, username, current_dir) #assuming for now that email is of format username@ucsd.edu


if __name__ == "__main__":
   main(sys.argv[1:])