project('Stable Distribution', 'cpp', default_options : ['cpp_std=c++11'],
         version : '0.9.7')

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

incdir=include_directories(['stable_distribution',
                            'Zolotarev'])

add_test_setup('gcc', env: 'LD_LIBRARY_PATH='+join_paths(meson.build_root(),'stable_distribution'))

subdir('doc')
subdir('stable_distribution')
subdir('basic_test')
subdir('fit_test')
subdir('quick_test')
subdir('FMStable_test')
subdir('dump_mpreal')
subdir('xcheck_to_file')
subdir('trace')
subdir('Zolotarev')
subdir('trace_zolotarev')
subdir('random_test')
subdir('duality_check')
