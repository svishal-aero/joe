project('joe', 'cpp')

sources_joe = [
  'common/src/Gp.cpp',
  'common/src/Logging.cpp',
  'common/src/MiscUtils.cpp',
  'common/src/MpiStuff.cpp',
  'common/src/MshFilter.cpp',
  'common/src/Param.cpp',
  'common/src/tc_vec3d.cpp',
  'common/src/Ugp.cpp',
  'common/src/UgpWithCv2.cpp',
  'common/src/UgpWithTools.cpp',
  'joe/src/JoeWithModels.cpp',
  'joe/src/JoeWithModelsAD.cpp',
  'joe/src/UgpWithCvCompFlow.cpp',
  'joe/src/UgpWithCvCompFlowAD.cpp',
  'joe/src/Scalars.cpp',
  'joe/src/ScalarsAD.cpp'
]

cxx = meson.get_compiler('cpp')

libPath = '/home/svishal/Git/joe/externals/petsc_install/lib'
colpack  = cxx.find_library('ColPack', dirs: [libPath])
adolc    = cxx.find_library('adolc', dirs: [libPath])
metis    = cxx.find_library('metis', dirs: [libPath])
parmetis = cxx.find_library('parmetis', dirs: [libPath])
petsc    = cxx.find_library('petsc', dirs: [libPath])

deps = [colpack, adolc, metis, parmetis, petsc]

flags = [
  '-I/home/svishal/Git/joe/externals/petsc_install/include',
  '-O3',
  '-w',
  '-fPIC',
  '-libverbs',
  '-fpermissive',
  '-std=c++11',
  '-Wno-format',
  '-Wl,--no-undefined',
  '-DMPI_OFFSET_IS_LONG_LONG_INT',
  '-DNO_ASSERT',
  '-DWITH_PARMETIS',
  '-DWITH_PETSC',
]

message('\nIgnoring CPP warnings\n')

joe = shared_library('Joe', sources_joe,
                     cpp_args: flags,
                     dependencies: deps,
                     install: true,
                     install_dir: 'lib',
                     install_rpath: libPath)
