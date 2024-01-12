import pyslha
import subprocess as sp
import numpy as np
import os
import threading
import shutil
import concurrent.futures

id = 'cMSSM'

# current working directory
SOURCE_DIR = os.getcwd()+'/'

# path to the micrOMEGAs binary
MO_PATH = SOURCE_DIR+'../../../external/micromegas_5.3.41/MSSM/scan_example'

# path to the SPheno binary
SPHENO_PATH = SOURCE_DIR+'../../../external/SPheno-4.0.5/bin/SPheno'

# create temporary folder
TEMP = 'temp_{}/'.format(id)
if not os.path.isdir(TEMP):
    os.mkdir(TEMP)

# default SPheno input file
DEFAULT_INPUT = 'cMSSM.in'
doc = pyslha.read(DEFAULT_INPUT,ignorenomass=True) # save the default input file in doc

# new SPheno input file
SPHENO_INPUT_FILE = SOURCE_DIR+TEMP+'SPheno_{}.in'

# SPheno output file
SPHENO_OUTPUT_FILE = SOURCE_DIR+TEMP+'LesHout_{}.spc'

# micrOMEGAs output file
MO_OUTPUT_FILE = SOURCE_DIR+TEMP+'MO_{}.out'


# spcification of the region of the parameter scan
N_0 = 200
N_12 = 200
m0_in = 3000
m12_in = 1400
m0 = np.linspace(2970,3040,N_0)
m12 = np.linspace(1380,1420,N_12)

# arrays to save results
Oh2MO,Oh2Tree,Oh2Full = -np.ones((N_0,N_12)),-np.ones((N_0,N_12)),-np.ones((N_0,N_12))
omegaCh_dmnlo= -np.ones((N_0,N_12))
mh,mst1,mneut1 = -np.ones((N_0,N_12)),-np.ones((N_0,N_12)),-np.ones((N_0,N_12))
pValMO,pValTree,pValFull = -np.ones((N_0,N_12)),-np.ones((N_0,N_12)),-np.ones((N_0,N_12))
csDDMO,csDDTree,csDDNLO = -np.ones((N_0,N_12,4)),-np.ones((N_0,N_12,4)),-np.ones((N_0,N_12,4))

# create folder for the results
RESULTS = os.path.abspath('results_{}'.format(id))
if not os.path.isdir(RESULTS):
    os.mkdir(RESULTS)



def calc(i,j):
    print(i,j)
    thread_id = threading.current_thread().name
    #print(thread_id)

    # modify soft-breaking parameters 
    doc.blocks['MINPAR'][1] = m0[j] # m0
    doc.blocks['MINPAR'][2] = m12[i] # m1/2

    # write new SPheno input
    pyslha.write(SPHENO_INPUT_FILE.format(thread_id),doc)

    # run SPheno to compute the spectrum
    output = sp.check_output([SPHENO_PATH,SPHENO_INPUT_FILE.format(thread_id), SPHENO_OUTPUT_FILE.format(thread_id)], stderr=sp.STDOUT, universal_newlines=True)
    if 'There has been a problem during the run.' in output or 'STOP' in output:
        return

    # save the higgs, stop and neutralino masses
    spectrum = pyslha.read(SPHENO_OUTPUT_FILE.format(thread_id),ignorenomass=True)
    mh[j,i] = spectrum.blocks['MASS'][25]
    mst1[j,i] = spectrum.blocks['MASS'][1000006]
    mneut1[j,i] = spectrum.blocks['MASS'][1000022]


    # run micrOMEGAs and get the relic density 
    call = sp.Popen([MO_PATH+' '+SPHENO_OUTPUT_FILE.format(thread_id)+' '+MO_OUTPUT_FILE.format(thread_id)], stdout=sp.PIPE, shell=True)
    stdout, stderr = call.communicate()
    print(stdout)
    if call.returncode != 0: 
        return

    results = pyslha.read(MO_OUTPUT_FILE.format(thread_id),ignorenomass=True)

    Oh2MO[j,i] = results.blocks['RESULTS'][1]
    Oh2Tree[j,i] = results.blocks['RESULTS'][2] # compute Oh2 with DM@NLO tree-level cross sections
    Oh2Full[j,i] = results.blocks['RESULTS'][3] # compute Oh2 with the full DM@NLO cross sections

    omegaCh_dmnlo[j,i] = results.blocks['RESULTS'][200]

    pValMO[j,i] = results.blocks['RESULTS'][41]
    pValTree[j,i] = results.blocks['RESULTS'][42]
    pValFull[j,i] = results.blocks['RESULTS'][43]

    for k in range(4):
        csDDMO[j,i,k] = results.blocks['RESULTS'][411+k]
        csDDTree[j,i,k] = results.blocks['RESULTS'][421+k]
        csDDNLO[j,i,k] = results.blocks['RESULTS'][431+k]


    with save_lock:

        os.chdir(RESULTS)

        # save values for future use
        np.save('m0.npy',m0)
        np.save('m12.npy',m12)

        np.save('mh.npy',mh)
        np.save('mstop1.npy',mst1)
        np.save('mneut1.npy',mneut1)

        np.save('Oh2MO.npy',Oh2MO)
        np.save('Oh2Tree.npy',Oh2Tree)
        np.save('Oh2Full.npy',Oh2Full)

        np.save('omegaCh_dmnlo.npy',omegaCh_dmnlo)

        np.save('pValMO.npy', pValMO)
        np.save('pValTree.npy',pValTree)
        np.save('pValFull.npy',pValFull)

        np.save('csDDMO.npy', csDDMO)
        np.save('csDDTree.npy',csDDTree)
        np.save('csDDNLO.npy',csDDNLO)

        os.chdir(SOURCE_DIR)


# Create a lock to ensure thread safety
save_lock = threading.Lock()

# Create and start the threads
max_threads = 36
with concurrent.futures.ThreadPoolExecutor(max_threads) as executor:
    for i in range(N_12): # loop over m0
        for j in range(N_0): # loop over m12
            executor.submit(calc,i,j)


# delete temporary folder
if os.path.isdir(TEMP):
    shutil.rmtree(TEMP) 


