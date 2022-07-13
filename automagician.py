#!/bin/env python

import inspect
import subprocess
import os
import time
import datetime
import gzip
from os.path import exists
from optparse import OptionParser

# This script does the following:

# 1. Record all jobs that is ever known. Containing date when calculation is recorded.
# 1.1 Service directory. Recursivly goes through every subdirectory and record it if it's not known.
# 1.2 At the time of registration, determine if the job is converged or not.

# 2. 

# defaults:
# Known jobs are stored in ~/known_jobs.dat
parser = OptionParser()
parser.add_option("-r","--register",action="store_true",dest="register",default=False,help="Process and add all calculations contained in the current directory and all its subdirectories to known_jobs.dat")
parser.add_option("-p","--process",action="store_true",dest="process",default=False,help="Check and re-process every job in unconverged_jobs.dat")
parser.add_option("-t","--test",action="store_true",dest="test",default=False,help="test the functions of this script without modifying original files or submitting jobs")
parser.add_option("-s","--silent",action="store_true",dest="silent",default=False,help="supress output")
parser.add_option("--rsc","--reset_converged",action="store_true",dest="reset_converged",default=False,help="move all converged jobs back into unconverged ones")
parser.add_option("--ac","--archive_converged",action="store_true",dest="archive_converged",default=False,help="Permanently log converged jobs")
parser.add_option("-l","--limit",action="store",dest="limit",type="int",default=99999,help="Stop this script if more than a number of jobs are in queue.")
parser.add_option("-m","--machine",action="store",dest="machine",type="int",default=0,help="Determines machine slot of submissions (Halifax, Stampeede, FRI), FRI used as default")
parser.parse_args()
print(parser.values)
#Chosen submission file
#0 - FRI, 1 - Halifax, 2 - Stampeede
def getsubfile():
    global subfile
    subfile = "frilab1.sub"
    if parser.values.machine == 1:
        subfile = "halifax.sub"
    elif parser.values.machine == 2:
        subfile = "skx.mpi.slurm"

  
# Look for jobs in the current directory and subdirectories to record.
def register():
  # os.walk() generates 3 tuples: the directory, subdirectories contained, files contained
  print("executed")
  target_files_list = ["POSCAR","POTCAR","INCAR","KPOINTS",subfile] 
  for stuff in os.walk(os.getcwd()):
    stop_limit(parser.values.limit)
    if not parser.values.silent :
      print("Registrator looking at ", stuff[0])
    # if a directory has those 5 files, it contains a VASP calculation
    has_calculation = True
    for target_file in target_files_list:
      if(target_file not in stuff[2]):
        #if not parser.values.silent :
        #  print("No calculations found")
        has_calculation = False
        break
      
    if(has_calculation):
      job_directory = stuff[0]
      known_jobs.seek(0)
      if(not job_directory in known_jobs.read()):
        record_job(job_directory)
      process_job(job_directory)

# add a job to the known_jobs.dat
def record_job(job_directory):
  known_jobs.seek(0)
  if( not ( job_directory in known_jobs.read())):
    known_jobs.write(job_directory + '\n') 

# This assumes that all converged calculations do not wrap up its last run
def determine_convergence(job_directory):
  # default
  # No CONTCAR and no ll_out is not converged
  if not ( exists(job_directory+"/CONTCAR") and exists(job_directory+"/ll_out") ):
    return False
  # use ll_out to determine convergence
  ll_out = open(job_directory+"/ll_out","r")
  ll_out.seek(0)
  os.chdir(job_directory)
  subprocess.call("vef.pl",stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
  if ("reached required accuracy - stopping structural energy minimisation" in ll_out.read()):
    ll_out.close()
    # ISIF = 3 requires that jobs converge after one step
    if special_instruction(job_directory) == 1 :
      if not parser.values.silent :
        print("This is a bulk relaxation job")
        return determine_box_convergence(job_directory)
    ll_out.close()
    return True
  else :
    ll_out.close()
    return False

def determine_box_convergence(job_directory):
  os.chdir(job_directory)
  fe = open("fe.dat",'r')
  fe_lines = fe.readlines()
  if (len(fe_lines) > 1 ):
    print("number of lines in fe.dat more than 1")
    fe.close()
    return False
  elif (len(fe_lines) == 1):
    if not parser.values.silent:
      print("box relaxation finished")
    fe.close()
    return True
  else:
    print("this calculation needs attention")
    fe.close()
    return False

# Determine if this job needs to be treated differently
def special_instruction(job_directory):
  # 1 --- ISIF = 3
  incar = open("INCAR",'r')
  for line in incar.readlines():
    if "ISIF" in line:
      line = line.strip("\n")
      if '3' in line.split('=')[1] :
        return 1
  return 0

def process_job(job_directory):
  stop_limit(parser.values.limit)
  if not parser.values.silent :
    print("processing job ",job_directory)
  is_converged = False
  is_running = running_check(job_directory)
  if is_running:
    if not parser.values.silent :
      print("job is running, do nothing")
  else:
    if(exists(job_directory+"/ll_out")):
      if check_error(job_directory) :
        record_error(job_directory)
        log_error(job_directory,home)
        fix_error(job_directory)
      if not exists(job_directory+"/CONTCAR"):
        qsub(job_directory)
      elif os.path.getsize(job_directory+"/CONTCAR") != 0:
        is_converged = determine_convergence(job_directory)
        if is_converged:
          print("job_converged")
          process_converged(job_directory)
        else:
          process_unconverged(job_directory)
      #qsub(job_directory)
    else :
      process_unconverged(job_directory)
      #qsub(job_directory)

# Check to see if a job has an error in it
def check_error(job_directory):
  os.chdir(job_directory)
  ll_out = open("ll_out",'r')
  if "I REFUSE TO CONTINUE WITH THIS SICK JOB" in ll_out.read():
    if not parser.values.silent:
      print("This job reported an error!")
    ll_out.close()
    return True
  else: 
    ll_out.close()
    return False

def record_error(job_directory):
  error_jobs.seek(0)
  if( not ( job_directory in error_jobs.read())):
    error_jobs.write(job_directory + '\n')

# generate a permanent error log  
def log_error(job_directory,home):
  error_log = open(home+"/error_log.dat","a+")
  for error_message in get_error_message(job_directory):
    error_log.write(str(datetime.datetime.now())+"  "+job_directory+"  "+error_message+"\n")



def get_error_message(job_directory):
  os.chdir(job_directory)
  ll_out = open('ll_out','r')
  messages=[]
  for line in ll_out:
    if ("ERROR" in line) or ("error" in line):
      messages.append(line)
  if len(messages) == 0:
    return "message not found!"  
  return messages

  

def fix_error(job_directory):
  os.chdir(job_directory)
  error_messages = get_error_message(job_directory)
  for error_message in error_messages:
    if parser.values.test:
      print(error_message)
      continue
    if "ZBRENT" in error_message :
      os.chdir(job_directory)
      if exists(job_directory+"/CONTCAR"):
        if os.path.getsize(job_directory+"/CONTCAR") != 0:
          wrap_up(job_directory)
      archive_ll_out()
      if not parser.values.silent:
        print("ll_out archived")
      qsub(job_directory)
      return None
    elif "number of potentials on File POTCAR incompatible with number" in error_message:
      os.chdir(job_directory)
      subprocess.call([home+'/kingRaychardsArsenal/sortpos.py'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
      subprocess.call([home+'/kingRaychardsArsenal/sogetsoftpbe.py'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
      archive_ll_out()
      qsub(job_directory)
      return None
  if not parser.values.silent:
    print("a fix was not attempted!")
    return None
  
def archive_ll_out():
  subprocess.call(['cat','ll_out','>','archive_stdout'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
  subprocess.call(['rm','ll_out'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)


def archive_converged(home):
  subprocess.call(['cat',home+'/converged_jobs.dat','>',home+'/archive_converged.dat'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
  subprocess.call(['rm',home+'/converged_jobs.dat'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

def reset_converged(home):
  subprocess.call(['cat',home+'/converged_jobs.dat','>',home+'/unconverged_jobs.dat'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
  subprocess.call(['rm',home+'/converged_jobs.dat'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

def process_converged(job_directory):
  # remove the job from unconverged_jobs.dat
  global unconverged_jobs
  unconverged_jobs.close()
  with unconverged_jobs_opener(home,'r') as f:
    old_lines = f.readlines()
  with unconverged_jobs_opener(home,'w') as f:
    for new_line in old_lines:
      if new_line.strip("\n") != job_directory:
        f.write(new_line)
  # add the job to converged_jobs.dat
  converged_jobs.seek(0)
  if not job_directory in converged_jobs.read():
    converged_jobs.write(job_directory+"\n")
  os.chdir(job_directory)
  subprocess.call("vef.pl",stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
  unconverged_jobs.close()
  unconverged_jobs = unconverged_jobs_opener(home,'a+')
  unconverged_jobs.seek(0)
  
def process_unconverged(job_directory):
    # First check if this is recorded as unconverged
    if not parser.values.silent:
      print("processing unconverged job at")
      print(job_directory)
    unconverged_jobs.seek(0)
    if not (job_directory in unconverged_jobs.read()):
      if not parser.values.silent:
        print("Adding this job to unconverged_jobs.dat")
      unconverged_jobs.write(job_directory+"\n")
    # is CONTCAR empty? empty file have size 0. if not empty, wrap up if job is not running
    if(not exists(job_directory+"/CONTCAR") or not exists(job_directory+"/OUTCAR")):
      qsub(job_directory)
    elif(os.path.getsize(job_directory+"/CONTCAR") != 0) :
      wrap_up(job_directory)
      qsub(job_directory)
    else:
      qsub(job_directory)
    # Then see if job is not in error
    # TODO
    # Then submit job if job is not already in the queue

        
def combine_XDAT_FE(job_directory):
  os.chdir(job_directory)
  cmbX = open("cmbXDATCAR","w")
  cmbF = open("cmbFE.dat","w")
  dirList = os.listdir("./")
  dirList = sorted(dirList)
  tlc = 0 # total line count
  # start from run0 and all the way up
  for dir in dirList :
   if os.path.isdir(dir) :
      # check XDATCAR.gz first
      if os.path.isfile(os.path.join(dir,"XDATCAR.gz")) :
        f = gzip.open(os.path.join(dir,"XDATCAR.gz"))
        for line in f.readlines() :
          try:
            line = line.decode("utf-8")
          except: 
            pass
          cmbX.write(line)
        f.close()
      # check for XDATCAR next
      if os.path.isfile(os.path.join(dir,"XDATCAR")) :
        f = open(os.path.join(dir,"XDATCAR"))
        for line in f.readlines() :
          try:
            line = line.decode("utf-8")
          except:
            pass
          cmbX.write(line)
        f.close()
      # now check for fe.dat
      if os.path.isfile(os.path.join(dir,"fe.dat")) :
        f = open(os.path.join(dir,"fe.dat"))
        for line in f.readlines() :
          cmbF.write(str(tlc) + "  " + line)
          tlc = tlc + 1
        f.close()
  # now check for XDATCAR at current directory
  if os.path.isfile("XDATCAR") :
    f = open("XDATCAR")
    for line in f.readlines() :
      cmbX.write(line + "\n")
    f.close()
  if os.path.isfile("fe.dat") :
    f = open("fe.dat")
    for line in f.readlines() :
      cmbF.write(str(tlc) + "  " + line)
      tlc = tlc + 1
    f.close()
  cmbX.close()
  cmbF.close()



  #finds out if a job is in queue
def running_check(job_directory):
  # first record all job status. 
  output = str(subprocess.check_output('qstat'))
  # if no job is running, output is b''
  #if( len(output) == 3):
  #  return False
  #else:
  job_list = output.split('\\n')
  for line in job_list[2:-1]:
    output = str(subprocess.check_output(['/opt/gridengine/bin/lx-amd64/qstat','-j',line.split()[0]]))
    qstat_j = output.split('\\n')
    for line in qstat_j:
      if("workdir" in line):
        if(line.split()[1] == job_directory):
          return True
  return False
  #is_running = True

def wrap_up(job_directory):
  if not parser.values.silent:
    print('wrapping up job')
  # first find name
  os.chdir(job_directory)
  directories = [f.path for f in os.scandir(job_directory) if f.is_dir()]
  runs = [ f for f in directories if "run" in f]
  # if no runs, wrap up into run0
  if( len(runs) == 0 ):
    if parser.values.test :
      print("this job would be wrapped up")
    else : 
      subprocess.call(['/usr/local/vtstscripts/vfin.pl','run0'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
  else:
    # find the largest run
    largest_number = 0
    for run in runs:
      number = int(run.partition('run')[2])
      if( number >= largest_number):
        largest_number = number + 1
    if not parser.values.test :
      subprocess.call(['/usr/local/vtstscripts/vfin.pl','run'+str(largest_number)],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
      combine_XDAT_FE(job_directory)
    else :
      print("this job would be wrapped into run "+str(largest_number))
  optimizer_review(job_directory)

# if CG isn't working, use DMD
def optimizer_review(job_directory):
  os.chdir(job_directory)
  INCAR = open("INCAR",'r')
  for INCAR_line in INCAR.readlines():
    if 'IBRION' in INCAR_line:
      INCAR_line = INCAR_line.strip("\n")
      if '2' in str(INCAR_line.split('=')[1]):
        if not os.path.exists("./cmbFE.dat"):
          combine_XDAT_FE(job_directory)
        cmbFE = open('cmbFE.dat','r')       
        cmbFE_lines = cmbFE.readlines()
        if len(cmbFE_lines) < 200:
          return None
        else :
          smallest = 0
          for i in range(1,10):
            #if float(cmbFE_lines('utf-8')[-i][4]) < smallest:
            try:
              cmbFE_lines[-i] = cmbFE_lines[-i].decode('utf-8')
            except:
              pass
            cmbFE_lines[-i] = cmbFE_lines[-i].split()
            if float(cmbFE_lines[-i][4]) < smallest:
              smallest = float(cmbFE_lines[-i][4])
            else:
              # CG may not be stable enough. change INCAR
              INCAR.seek(0)
              INCAR_lines = INCAR.readlines()
              INCAR.close()
              if parser.values.test:
                print("INCAR will be modified")
                INCAR = open('testINCAR','w')
                for INCAR_line in INCAR_lines:
                  if 'IBRION' not in INCAR_line:
                    INCAR.write(INCAR_line)
                  else:
                    INCAR.write("IBRION = 3 \n")
                INCAR.close()
                return None
              else:
                INCAR = open('INCAR','w')
                for INCAR_line in INCAR_lines:
                  if 'IBRION' not in INCAR_line:
                    INCAR.write(INCAR_line)
                  else:
                    INCAR.write("IBRION = 3 \n")
                INCAR.close()
                return None
              
      else:
        return None
      
  

def qsub(job_directory):
  if parser.values.test:
    print("qsub is called by ", inspect.stack()[1].function)
    print("a job would have been submitted here")
  else:
    os.chdir(job_directory)
    update_job_name()
    print("submitting job")
    if parser.values.machine == 0 or parser.values.machine == 1:
        subprocess.run(["qsub",subfile])
    elif parser.values.machine == 2:
        subprocess.run(["sbatch",subfile])
    else:
        print("Unknown submission method")

def update_job_name():
  script = open(subfile,"r")
  script_lines = script.readlines()
  script.close()
  script = open(subfile,"w")
  for line in script_lines:
    print(line)
    if "-N" in line:
      script.write("#$ -N "+"path_"+ os.getcwd().replace("/","_")+"\n")
    else:
      script.write(line)

def stop_limit(limit):
  output = str(subprocess.check_output('qstat'))
  job_list = output.split('\\n')
  if len(job_list) >= limit :
    print("Hit job limit of ", parser.values.limit)
    print("Exiting.")
    exit()

def update_known():
  unconverged_jobs.seek(0)
  job_list = [ entry.strip("\n") for entry in unconverged_jobs]
  for job_directory in job_list:
    process_job(job_directory)
    
def known_jobs_opener(home,mode):
  if parser.values.test:
    subprocess.call(['cp',home+'/known_jobs.dat',home+'/test_known_jobs.dat'])
    known_jobs = open(home+'/test_known_jobs.dat',mode)
  else :
    known_jobs = open(home+'/known_jobs.dat',mode)
  return known_jobs

def unconverged_jobs_opener(home,mode):
  if parser.values.test :
    subprocess.call(['cp',home+'/unconverged_jobs.dat',home+'/test_unconverged_jobs.dat'])
    unconverged_jobs = open(home+'/test_unconverged_jobs.dat',mode)
  else :
    unconverged_jobs = open(home+'/unconverged_jobs.dat',mode)
  return unconverged_jobs

def converged_jobs_opener(home,mode):
  if parser.values.test :
    subprocess.call(['cp',home+'/converged_jobs.dat',home+'/test_converged_jobs.dat'])
    converged_jobs = open(home+'/test_converged_jobs.dat',mode)
  else :
    converged_jobs = open(home+'/converged_jobs.dat',mode)
  return converged_jobs

def error_jobs_opener(home,mode):
  if parser.values.test :
    subprocess.call(['cp',home+'/error_jobs.dat',home+'/test_error_jobs.dat'])
    error_jobs = open(home+'/test_error_jobs.dat',mode)
  else:
    error_jobs = open(home+'/error_jobs.dat',mode)
  return error_jobs

getsubfile()
home = os.environ['HOME']
known_jobs = known_jobs_opener(home,'a+')
known_jobs.seek(0)
unconverged_jobs = unconverged_jobs_opener(home,'a+')
unconverged_jobs.seek(0)
converged_jobs = converged_jobs_opener(home,'a+')
converged_jobs.seek(0)
error_jobs = error_jobs_opener(home,'a+')
error_jobs.seek(0)

if parser.values.reset_converged:
  reset_converged(home)
if parser.values.archive_converged:
  archive_converged(home)

if parser.values.register:
  print("Registerng all jobs in the current directory")
  register()
if parser.values.process:
  print("Processing all jobs in unconverged_jobs.dat")
  update_known()


