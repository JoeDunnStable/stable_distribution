This directory contains .pc files that are needed by meson and that
were missing in one or more of platforms.  They should be edited and 
placed in the pkgconfig directory where meson can find them.

The mpreal.h file is slightly modified from the original to allow 
use with boost by commenting out the define 

MPREAL_HAVE_DYNAMIC_STD_NUMERIC_LIMITS 

and setting digits = 96, digits10 = 28, and max_digits10=29.
