#Author:Chayan
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC
import math, re
import argparse
import ftplib
import socket
import random
import time
from random import randint
import colorsys
import os, sys, os.path, math
import gzip
import getopt
from collections import OrderedDict
import subprocess
import glob
import textwrap

usage= ''' Description:  HmmSearch per Protein Accession '''
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-a", "--accession", help=" Protein Accession ")
parser.add_argument("-f", "--fasta", help=" Protein fasta ")
parser.add_argument("-d", "--localHmmDirectory", help=" Path for Local Files, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-e", "--ethreshold", help=" E value threshold. Default = 1e-10 ")
parser.add_argument("-u", "--user_email", required=True, action="append", metavar="RECIPIENT",default=[], dest="recipients", help=" User Email Address (at least one required) ")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.0')
parser.add_argument("-k", "--keep", action="store_true", help=" Use this option to keep the intermediate files. ")
args = parser.parse_args()
parser.parse_args()

if args.ethreshold:
	evthresh=args.ethreshold
else:
	evthresh="1e-10"

ncbi_time=0.5
timeout = 10
socket.setdefaulttimeout(timeout)

Entrez.email = args.recipients[0] #User email
Entrez.max_tries = 5
Entrez.sleep_between_tries = 60


if args.localHmmDirectory:
	if os.path.isdir(args.localHmmDirectory):
		if args.localHmmDirectory[-1]=='/':
			localH=args.localHmmDirectory
			print('Local HMMs Data path : ', localH, '\n')
		else:
			localH=args.localHmmDirectory+'/'
			print('Local HMMs Data path : ', localH, '\n')
	else:
		print('No directory Found as : '+ args.localHmmDirectory)
		sys.exit()
else:
	localH='./'

#print(localH)

def seq_from_wp(accession_nr):
	"""
	:param accession_nr: NCBI protein accession
	:return: Protein Sequence
	"""
	try:
		time.sleep(ncbi_time)
		handle = Entrez.efetch(db="protein", id=accession_nr, rettype="gbwithparts", retmode="text")
	except Exception as e:
		print(str(e), ", error in entrez-fetch protein accession, {}, not found in database. \n" "Continuing with the next protein in the list. \nError in function: {}".format(accession_nr, seq_from_wp.__name__))
		return False

	record = SeqIO.read(handle, "genbank")
	handle.close()
	return record.format("fasta")

dirName=''
if args.accession:
	dirName=args.accession
else:
	if '.fasta'==args.fasta[-6:]:
		dirName=args.fasta[:-6]
	if '.faa'==args.fasta[-4:]:
		dirName=args.fasta[:-4]
	else:
		print('Sequence file should have extension like ".fasta" or ".faa"')
		sys.exit()


HmmOutDir=dirName+'_outDir'
if not os.path.exists(dirName+'_outDir'):
	os.makedirs(dirName+'_outDir')

if args.accession:
	with open (args.accession+'.faa','w') as fileOut:
		print(seq_from_wp(args.accession), file=fileOut)

'''
if not args.accession:
	fasIn=open(args.fasta, 'r').read().rstrip().split('>')
	with open(dirName+'.faa', 'w') as fileOut:
		for items in fasIn:
			if items!='':
				fastaID='>'+items.split('\n')[0]
				sequence=textwrap.fill(items.split('\n')[1],80)
				print(fastaID,sequence,sep='\n', file=fileOut)
'''

hmmList=[]
for hmms in (glob.glob(localH+"*.hmm")):
	hmmList.append(hmms)

#print(hmmList)

tabList=[]

for items in hmmList:
	tabFile=HmmOutDir+'/'+items.split('/')[-1][:-4]+'_tab'
	tabList.append(tabFile)
	alnFile=HmmOutDir+'/'+items.split('/')[-1][:-4]+'_aln'
	outFile=HmmOutDir+'/'+items.split('/')[-1][:-4]+'_out'
	if args.accession:
		command="hmmsearch --domtblout %s -A %s --incE %s --incdomE %s -E %s %s %s > %s" %(tabFile, alnFile, evthresh, evthresh, evthresh, items, './'+args.accession+'.faa', outFile)
	else:
		command="hmmsearch --domtblout %s -A %s --incE %s --incdomE %s -E %s %s %s > %s" %(tabFile, alnFile, evthresh, evthresh, evthresh, items, './'+args.fasta, outFile)
	#print(command)
	os.system(command)


ieList=[]
ieval=[]
for files in tabList:
	with open(files, 'r') as tabIn:
		for line in tabIn:
			Line=line.rstrip()
			LineSplit=re.sub('\s+', '\t', Line).split('\t')
			if Line[0]!='#':
				ieList.append(LineSplit)
				ieval.append(float(LineSplit[12]))
				#print(LineSplit[0], LineSplit[3], LineSplit[12])

if len(ieList)>0:
	for items in ieList:
		#print(items)
		if float(items[12])==min(ieval):
			print('# Query Protein Accession: '+items[0]+'\n'+'-- HMM match: '+items[3]+'\n'+'-- Independent Evalue: '+ items[12]+'\n')
else:
	print('No hits found with the current parameters')

if not args.keep:
	import shutil
	if os.path.exists(dirName+'_outDir'):
		shutil.rmtree(dirName+'_outDir')
