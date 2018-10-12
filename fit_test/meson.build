

fit_test = executable('fit_test', 
                        'fit_test.cpp',
                        dependencies : stable_distribution_dep) 


my_timeout = 4000

test('fit_test', fit_test, args: 'q_mle', timeout: my_timeout)

