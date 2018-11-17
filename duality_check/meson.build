

duality_check = executable('duality_check', 
                        'duality_check.cpp',
                        dependencies : stable_distribution_dep) 


my_timeout = 10000

test('duality_check_double', duality_check, args : 'double', timeout: my_timeout)
test('duality_check_mpreal', duality_check, args : 'mpreal', timeout: my_timeout)
test('duality_check_mpfr_float', duality_check, args : 'mpfr_float', timeout: my_timeout)
test('duality_check_cpp_bin_float', duality_check, args : 'cpp_bin_float', timeout: my_timeout)

