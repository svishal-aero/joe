import os, sys
from shutil import which
from subprocess import call as subp_call

def call(cmd):
    subp_call(cmd, shell=True)

def getOptions():
    options = []
    options.extend(sys.argv[1:])
    if len(options)==0:
        options.extend(['download', 'configure', 'make', 'install'])
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
    homeDir   = os.environ['HOME']
    workDir   = os.getcwd()
    thisDir   = os.path.dirname(os.path.abspath(__file__))
    petscDir  = os.path.join(thisDir, 'petsc')
    prefixDir = os.path.join(thisDir, 'petsc_install')

    call('mkdir -p '+os.path.basename(prefixDir))
    call('rm -rf '+os.path.join(prefixDir,'*'))

    os.chdir(thisDir)
    
    if 'download' in options:
        if not os.path.exists(petscDir):
            call('git clone -b release https://gitlab.com/petsc/petsc.git petsc')

    os.chdir(petscDir)
    
    if 'configure' in options:
        config_cmd  = './configure'
        config_cmd += ' --prefix='+prefixDir
        config_cmd += ' --with-mpi-dir='+mpiDir
        config_cmd += ' --download-colpack=1'
        config_cmd += ' --download-adolc=1'
        config_cmd += ' --download-f2cblaslapack=1'
        config_cmd += ' --download-metis=1'
        config_cmd += ' --download-parmetis=1'
        print('Running command: '+config_cmd)
        call(config_cmd)

    if 'make' in options:   
        call('make')

    if 'install' in options:
        call('make install')
    
    os.chdir(workDir)
