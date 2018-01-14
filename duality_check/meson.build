

duality_check = executable('duality_check', 
                        'duality_check.cpp',
                        dependencies : stable_distribution_dep) 


my_timeout = 100

test('duality_check_double', duality_check, args : 'double', timeout: my_timeout)

