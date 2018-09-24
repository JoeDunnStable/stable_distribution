project('Stable Distribution', 'cpp', default_options : ['cpp_std=c++11'],
         version : '1.0.2')

compiler=meson.get_compiler('cpp').get_id()

if compiler == 'gcc' or compiler == 'clang'
  add_global_arguments(
      [ '-fvisibility=hidden' ,
	       '-Wno-non-virtual-dtor' ,
	       '-Wno-sign-compare' ,
	       '-Wno-uninitialized',
               '-Wno-return-type' ,
               '-Wno-parentheses'
      ] ,
      language : 'cpp' )
elif compiler == 'msvc'
  add_global_arguments(
      ['/bigobj',
       '/wd4146',
       '/wd4244' ,
       '/wd4251' ,
       '/wd4715'
      ] ,
      language : 'cpp' )
  add_global_link_arguments(
       '/STACK:16000000',
      language : 'cpp' )
endif

#Eigen3 headers
eigendep = dependency('eigen3', version : '>3.0', required : true)

#Boost header and some libraries.
#I had to define BOOST_ROOT to get this to work.
boostdep = dependency('boost', modules : ['system', 'filesystem', 'timer'],  required : true)

# The Gnu multiprecision library
mpfrdep = dependency('mpfr', required : true)

# Pavel Holoborodko's C++ wrapper for gnu's mpfr
mpfrcppdep = dependency('mpfrcpp', required : true)

alldep = [eigendep, boostdep, mpfrdep, mpfrcppdep]

if host_machine.system() == 'linux'
  alldep = [alldep, meson.get_compiler('cpp').find_library('pthread')]
endif


stable_config = configuration_data()
stable_config.set('PACKAGE_NAME', meson.project_name())
stable_config.set('PACKAGE_VERSION', meson.project_version())
stable_config.set('PACKAGE_COMPILER', meson.get_compiler('cpp').get_id())
stable_config.set('PACKAGE_OS',host_machine.system())

stable_config_h = configure_file(input : 'stable_config.h.in',
				  output : 'stable_config.h',
				  configuration : stable_config,
				  install : false)

incdir=include_directories(['.',
                            'stable_distribution',
                            'Zolotarev'])

add_test_setup('gcc', env: 'DYLD_LIBRARY_PATH='+join_paths(meson.build_root(),'stable_distribution'))

subdir('doc')
subdir('stable_distribution')
subdir('basic_test')
subdir('fit_test')
subdir('quick_test')
subdir('FMStable_test')
subdir('stable_dump')
subdir('xcheck_to_file')
subdir('trace')
subdir('Zolotarev')
subdir('trace_zolotarev')
subdir('random_test')
subdir('duality_check')
