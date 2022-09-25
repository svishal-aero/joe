import os, sys
from shutil import which
from subprocess import call as subp_call

def call(cmd):
    subp_call(cmd, shell=True)

def getOptions():
    options = []
    options.extend(sys.argv[1:])
    if len(options)==0:
        options.extend(['petsc', 'metis', 'parmetis'])
    return options

def getMpiDir():
    mpirun_path = which('mpirun') 
    if mpirun_path is not None:
        mpiDir = os.path.dirname(os.path.dirname(os.path.realpath(mpirun_path)))
        print('Using MPI_DIR='+mpiDir)
    else:
        assert False, 'Please install an MPI implementation (e.g. openmpi, mpich)'
    return mpiDir

if __name__=='__main__':

    options   = getOptions()

    mpiDir    = getMpiDir()
    workDir   = os.getcwd()
    thisDir   = os.path.dirname(os.path.abspath(__file__))
    prefixDir = os.path.join(thisDir, 'opt')

    call('mkdir -p '+os.path.basename(prefixDir))
    call('rm -rf '+os.path.join(prefixDir,'*'))

    if 'petsc' in options:
        os.chdir(thisDir)
        call('tar -xzvf petsc-3.17.tar.gz')
        os.chdir('petsc-3.17.4')
        call('./configure --prefix='+prefixDir+' --with-mpi-dir='+mpiDir+' --download-colpack=1 --download-adolc=1 --download-f2cblaslapack=1 --with-fc=0')
        call('make')
        call('make install')
        os.chdir(workDir)

    if 'metis' in options:
        os.chdir(thisDir)
        call('tar -xzvf metis-5.1.0.tar.gz')
        os.chdir('metis-5.1.0')
        call('make config shared=1 prefix='+prefixDir)
        call('make -j')
        call('make install')
        os.chdir(workDir)

    if 'parmetis' in options:
        os.chdir(thisDir)
        call('tar -xzvf parmetis-4.0.3.tar.gz')
        os.chdir('parmetis-4.0.3')
        call('make config shared=1 prefix='+prefixDir)
        call('make -j')
        call('make install')
        os.chdir(workDir)
