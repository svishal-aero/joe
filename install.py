# There are three options
# "python install.py" will perform any required installation steps to build joe libraries
# "python install.py clean" will remove all build files but retain external dependencies
# "python install.py purge" will remove all build files and external dependencies (petsc, etc.)

import os, sys
from subprocess import call as subp_call

def call(cmd):
    subp_call(cmd, shell=True)


def getCC():
    #return which('mpicxx')
    return 'mpicxx'


def getCFLAGS(petsc_include):
    cflags  = ''
    cflags += ' -I'+petsc_include
    cflags += ' -O3'
    cflags += ' -fPIC'
    cflags += ' -libverbs'
    cflags += ' -fpermissive'
    cflags += ' -std=c++11'
    cflags += ' -Wno-format'
    cflags += ' -Wl,--no-undefined'
    cflags += ' -DMPI_OFFSET_IS_LONG_LONG_INT'
    cflags += ' -DWITH_PARMETIS'
    cflags += ' -DNO_ASSERT'
    cflags += ' -DWITH_PETSC'
    return cflags


def getCLIBS(petsc_lib):
    clibs  = ''
    clibs += ' -L'+petsc_lib
    clibs += ' -Wl,-rpath='+petsc_lib
    clibs += ' -lColPack'
    clibs += ' -ladolc'
    clibs += ' -lmetis'
    clibs += ' -lparmetis'
    clibs += ' -lpetsc'
    return clibs


def write_joecxx(joeHome, petsc_include, petsc_lib):

    call('rm -f joecxx')

    includes  = ''
    includes += ' -I'+os.path.join(joeHome, 'joe', 'include')
    includes += ' -I'+os.path.join(joeHome, 'common', 'include')
    includes += getCFLAGS(petsc_include)

    libs  = ''
    libs += getCLIBS(petsc_lib)
    libs += ' -L'+os.path.join(joeHome, 'common', 'lib')
    libs += ' -Wl,-rpath='+os.path.join(joeHome, 'common', 'lib')
    libs += ' -lJoeCommon'
    libs += ' -L'+os.path.join(joeHome, 'joe', 'lib')
    libs += ' -Wl,-rpath='+os.path.join(joeHome, 'joe', 'lib')
    libs += ' -lJoe'

    with open('joecxx', 'w') as f:
        f.write('#!'+sys.executable+'\n')
        f.write('import sys\n')
        f.write('from subprocess import call\n')
        f.write('call("mpicxx'
               +includes
               +' "+(" ").join(sys.argv[1:])+"'
               +libs
               +'", shell=True)\n')

    call('chmod +x joecxx')


def install():

    homeDir = os.environ['HOME']
    workDir = os.getcwd()
    joeHome = os.path.dirname(os.path.abspath(__file__))

    os.chdir(joeHome)

    petsc_dir     = os.path.join(joeHome, 'externals', 'petsc_install')
    petsc_lib     = os.path.join(petsc_dir, 'lib')
    petsc_include = os.path.join(petsc_dir, 'include')

    if not os.path.exists(petsc_dir):
        print()
        print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print(r' JOE INSTALLER: Installing external packages')
        print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        os.chdir('externals')
        call('python install.py')
        os.chdir('..')

    with open('Makefile.in', 'w') as f:
        print()
        print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print(r' JOE INSTALLER: Writing "Makefile.in"')
        print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        f.write('CC =%s' % getCC())
        f.write('\n\nCFLAGS =%s' % getCFLAGS(petsc_include))
        f.write('\n\nCLIBS =%s' % getCLIBS(petsc_lib))

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Calling "make"')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    call('make -j')

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Creating "joecxx"')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    write_joecxx(joeHome, petsc_include, petsc_lib)
    os.chdir(os.path.join(homeDir, 'bin'))
    call('rm -f joecxx')
    call('ln -s '+os.path.join(joeHome,'joecxx')+' .')

    os.chdir(workDir)


def clean():

    workDir = os.getcwd()
    joeHome = os.path.dirname(os.path.abspath(__file__))

    petsc_dir     = os.path.join(joeHome, 'externals', 'petsc_install')
    petsc_lib     = os.path.join(petsc_dir, 'lib')
    petsc_include = os.path.join(petsc_dir, 'include')

    if not os.path.exists('Makefile.in'):
        with open('Makefile.in', 'w') as f:
            print()
            print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            print(r' JOE INSTALLER: Writing "Makefile.in"')
            print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            f.write('CC =%s' % getCC())
            f.write('\n\nCFLAGS =%s' % getCFLAGS(petsc_include))
            f.write('\n\nCLIBS =%s' % getCLIBS(petsc_lib))

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Calling "make clean"')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    os.chdir(joeHome)
    call('make clean')
    os.chdir(workDir)


def purge():

    homeDir = os.environ['HOME']
    workDir = os.getcwd()
    joeHome = os.path.dirname(os.path.abspath(__file__))

    os.chdir(joeHome)

    clean()

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Removing "Makefile.in"')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    call('rm -f Makefile.in')

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Removing external packages')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    os.chdir('externals')
    call('rm -rf petsc*')
    os.chdir('..')

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Removing "joecxx"')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    os.chdir(os.path.join(homeDir, 'bin'))
    call('rm -f joecxx')
    os.chdir(joeHome)
    call('rm -f joecxx')

    os.chdir(workDir)


if __name__=='__main__':

    if len(sys.argv)==1:       install()
    elif sys.argv[1]=='purge': purge()
    elif sys.argv[1]=='clean': clean()
    else: print('Error: Unrecognized option '+sys.argv[1])
