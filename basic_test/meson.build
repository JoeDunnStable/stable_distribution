

basic_test = executable('basic_test', 
                        ['basic_test.cpp', 'kolmogorov_cdf.cpp', 'kolmogorov.h'],
                        dependencies : stable_distribution_dep) 


double_timeout = 100

test('cdf_double', basic_test, args: 'cdf_double', timeout: double_timeout)
test('pdf_double', basic_test, args: 'pdf_double', timeout: double_timeout)
test('ddx_pdf_double', basic_test, args: 'ddx_pdf_double', timeout: double_timeout)
test('quantile_double', basic_test, args: 'quantile_double', timeout: double_timeout)
test('mode_double', basic_test, args: 'mode_double', timeout: double_timeout)
test('random_double', basic_test, args: 'random_double', timeout: double_timeout)

mpreal_timeout = 4000

test('cdf_mpreal', basic_test, args: 'cdf_mpreal', timeout: mpreal_timeout)
test('pdf_mpreal', basic_test, args: 'pdf_mpreal', timeout: mpreal_timeout)
test('ddx_pdf_mpreal', basic_test, args: 'ddx_pdf_mpreal', timeout: mpreal_timeout)
test('quantile_mpreal', basic_test, args: 'quantile_mpreal', timeout: mpreal_timeout)
test('mode_mpreal', basic_test, args: 'mode_mpreal', timeout: mpreal_timeout)
test('random_mpreal', basic_test, args: 'random_mpreal', timeout: mpreal_timeout)

mpfr_float_timeout = 5000

test('cdf_mpfr_float', basic_test, args: 'cdf_mpfr_float', timeout: mpfr_float_timeout)
test('pdf_mpfr_float', basic_test, args: 'pdf_mpfr_float', timeout: mpfr_float_timeout)
test('ddx_pdf_mpfr_float', basic_test, args: 'ddx_pdf_mpfr_float', timeout: mpfr_float_timeout)
test('quantile_mpfr_float', basic_test, args: 'quantile_mpfr_float', timeout: mpfr_float_timeout)
test('mode_mpfr_float', basic_test, args: 'mode_mpfr_float', timeout: mpfr_float_timeout)
test('random_mpfr_float', basic_test, args: 'random_mpfr_float', timeout: mpfr_float_timeout)

cpp_bin_float_timeout = 6000

test('cdf_cpp_bin_float', basic_test, args: 'cdf_cpp_bin_float', timeout: cpp_bin_float_timeout)
test('pdf_cpp_bin_float', basic_test, args: 'pdf_cpp_bin_float', timeout: cpp_bin_float_timeout)
test('ddx_pdf_cpp_bin_float', basic_test, args: 'ddx_pdf_cpp_bin_float', timeout: cpp_bin_float_timeout)
test('quantile_cpp_bin_float', basic_test, args: 'quantile_cpp_bin_float', timeout: cpp_bin_float_timeout)
test('mode_cpp_bin_float', basic_test, args: 'mode_cpp_bin_float', timeout: cpp_bin_float_timeout)
test('random_cpp_bin_float', basic_test, args: 'random_cpp_bin_float', timeout: cpp_bin_float_timeout)

