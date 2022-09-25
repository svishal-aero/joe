# There are three options
# "python install.py" will perform any required installation steps to build joe libraries
# "python install.py clean" will remove all build files but retain external dependencies
# "python install.py purge" will remove all build files and external dependencies (petsc, etc.)

import os, sys
from shutil import which
from subprocess import call as subp_call

def call(cmd):
    subp_call(cmd, shell=True)


def getCFLAGS():
    cflags  = ''
    cflags += ' -O3'
    cflags += ' -fpermissive'
    cflags += ' -std=c++11'
    cflags += ' -Wno-format'
    cflags += ' -Wl,--no-undefined'
    cflags += ' -DMPI_OFFSET_IS_LONG_LONG_INT'
    cflags += ' -DWITH_PARMETIS'
    cflags += ' -DNO_ASSERT'
    cflags += ' -DWITH_PETSC'
    return cflags


def write_joecxx(joeHome):

    call('rm -f joecxx')

    includes  = ''
    includes += ' -I'+os.path.join(joeHome, 'joe', 'include')
    includes += getCFLAGS()

    libs  = ''
    libs += ' -L'+os.path.join(joeHome, 'install', 'lib')
    libs += ' -Wl,-rpath='+os.path.join(joeHome, 'install', 'lib')
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

    opt_dir = os.path.join(joeHome, 'externals', 'opt')

    if not os.path.exists(opt_dir):
        print()
        print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print(r' JOE INSTALLER: Installing external packages')
        print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        os.chdir('externals')
        call('python install.py')
        os.chdir('..')

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Writing meson.build file')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    if not os.path.exists('meson.build'):
        with open('meson.build', 'w') as f:
            f.write("project('joe', 'cpp')\n\n")

            f.write("sources_joe = [\n")
            f.write("  'common/src/Gp.cpp',\n")
            f.write("  'common/src/Logging.cpp',\n")
            f.write("  'common/src/MiscUtils.cpp',\n")
            f.write("  'common/src/MpiStuff.cpp',\n")
            f.write("  'common/src/MshFilter.cpp',\n")
            f.write("  'common/src/Param.cpp',\n")
            f.write("  'common/src/tc_vec3d.cpp',\n")
            f.write("  'common/src/Ugp.cpp',\n")
            f.write("  'common/src/UgpWithCv2.cpp',\n")
            f.write("  'common/src/UgpWithTools.cpp',\n")
            f.write("  'joe/src/JoeWithModels.cpp',\n")
            f.write("  'joe/src/JoeWithModelsAD.cpp',\n")
            f.write("  'joe/src/UgpWithCvCompFlow.cpp',\n")
            f.write("  'joe/src/UgpWithCvCompFlowAD.cpp',\n")
            f.write("  'joe/src/Scalars.cpp',\n")
            f.write("  'joe/src/ScalarsAD.cpp'\n")
            f.write("]\n\n")

            f.write("cxx = meson.get_compiler('cpp')\n\n")

            f.write("libPath = '"+opt_dir+"/lib'\n")

            f.write("colpack  = cxx.find_library('ColPack',  dirs: [libPath])\n")
            f.write("adolc    = cxx.find_library('adolc',    dirs: [libPath])\n")
            f.write("metis    = cxx.find_library('metis',    dirs: [libPath])\n")
            f.write("parmetis = cxx.find_library('parmetis', dirs: [libPath])\n")
            f.write("petsc    = cxx.find_library('petsc',    dirs: [libPath])\n\n")

            f.write("deps = [colpack, adolc, parmetis, petsc]\n\n")

            f.write("flags = [\n")
            f.write("  '-I"+opt_dir+"/include',\n")
            f.write("  '-O3',\n")
            f.write("  '-w',\n")
            f.write("  '-fPIC',\n")
            f.write("  '-libverbs',\n")
            f.write("  '-fpermissive',\n")
            f.write("  '-std=c++11',\n")
            f.write("  '-Wno-format',\n")
            f.write("  '-Wl,--no-undefined',\n")
            f.write("  '-DMPI_OFFSET_IS_LONG_LONG_INT',\n")
            f.write("  '-DNO_ASSERT',\n")
            f.write("  '-DWITH_PARMETIS',\n")
            f.write("  '-DWITH_PETSC',\n")
            f.write("]\n\n")

            f.write("joe = shared_library('Joe', sources_joe,\n")
            f.write("                     cpp_args: flags,\n")
            f.write("                     dependencies: deps,\n")
            f.write("                     install: true,\n")
            f.write("                     install_dir: 'lib',\n")
            f.write("                     install_rpath: libPath)\n")

    if not os.path.exists('build'):
        call('CXX=mpicxx meson setup build --prefix='+joeHome+'/install')
    call('meson compile -C build')
    call('meson install -C build --only-changed')

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Creating "joecxx"')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    write_joecxx(joeHome)
    os.chdir(os.path.join(homeDir, 'bin'))
    call('rm -f joecxx')
    call('ln -s '+os.path.join(joeHome, 'joecxx')+' .')

    os.chdir(workDir)


def clean():

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Removing "build" and "install" directories')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    call('rm -rf build install')


def purge():

    homeDir = os.environ['HOME']
    workDir = os.getcwd()
    joeHome = os.path.dirname(os.path.abspath(__file__))

    os.chdir(joeHome)

    clean()

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Removing meson.build file')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    call('rm -f meson.build')

    print()
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print(r' JOE INSTALLER: Removing external packages')
    print(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    os.chdir('externals')
    call('rm -rf parmetis-4.0.3 metis-5.1.0 petsc-3.17.4')
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

    assert(which('mpicxx') is not None)
    assert(which('meson') is not None)
    if len(sys.argv)==1:       install()
    elif sys.argv[1]=='purge': purge()
    elif sys.argv[1]=='clean': clean()
    else: print('Error: Unrecognized option '+sys.argv[1])
