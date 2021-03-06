import os
import time
from subprocess import *

nodes = [512] #, 1024, 2048, 4096, 8192, 16384, 32768]

def runcmd (node):

	script = './runbgq.sh ' + str(node) 
	cmd = 'qsub -A Performance -t 01:00:00 -n '+str(node)+' --mode script '+script
	print 'Executing ' + cmd
	jobid = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
	print 'Jobid : ' + jobid
	
	while True:
		cmd = 'qstat ' + jobid.strip() + ' | grep preeti | awk \'{print $1}\''
		jobrun = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
		if jobrun == '':
			break
		time.sleep(30)

	return jobid.strip()

for iter in range (2, 3):
 for node in nodes:
		print '\nStarting on ' + str(node) + ' nodes' #+ str(rank) + ' ranks per node'
		jobid = runcmd(node)
		filename = 'rt_'+str(node)+'_'+str(iter)
		print filename + ' ' + jobid
		cmd = 'mv ' + jobid.strip() + '.output ' + filename + '.output'
		Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
		cmd = 'mv ' + jobid.strip() + '.error ' + filename + '.error'
		Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
		cmd = 'mv ' + jobid.strip() + '.cobaltlog ' + filename + '.cobaltlog'
		Popen(cmd, shell=True, stdout=PIPE).communicate()[0]

		cmd = 'ls -tr it_N*'
		output_files = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
		opfname = output_files.split()
		print opfname

		cmd = 'ls -tr mpi_profile*.0'
		mpifnames = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
		mpifname = mpifnames.split()
		print mpifname

		for i in range(len(opfname)):
			print i
			cmd = 'mv ' + mpifname[i].strip() + ' ' + opfname[i].strip() + '-' + mpifname[i].strip()
			print cmd
			Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
			cmd = 'mv ' + opfname[i].strip() + ' o-' + opfname[i].strip()
			Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
			print cmd

