
import numpy.distutils.intelccompiler
from setuptools import setup, Extension, find_packages
import os


#-D_GLIBCXX_USE_CXX11_ABI=0: because GCC 5 issue with dual ABI: https://gcc.gnu.org/onlinedocs/gcc-5.2.0/libstdc++/manual/manual/using_dual_abi.html

def extcpp(name,
        sources=[],
        include_dirs=['utils/ext/', os.environ['ANACONDA_ENV']+'/include/python2.7', os.environ['ANACONDA_ENV']+'/include', '/usr/local/include'],
        library_dirs=[os.environ['ANACONDA_ENV']+'/lib'],
        libraries=[],
        extra_compile_args=['-Wall -g -shared -O2 -fPIC -fopenmp -Wl,--export-dynamic -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0'],
        extra_link_args=['-lpython2.7 -lboost_python -lboost_iostreams -lpthread -larmadillo -lmsgpack -fopenmp']):
    return Extension(name,
                     include_dirs=include_dirs,
                     library_dirs=library_dirs,
                     libraries=libraries,
                     sources=sources,
                     extra_compile_args=extra_compile_args,
                     extra_link_args=extra_link_args)

def extc(name, sources=[], include_dirs=[], library_dirs=[], libraries=[], extra_compile_args=['-g -fopenmp -std=c99'], extra_link_args=['-g -fopenmp']):
    return Extension(name, include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries, sources=sources, extra_compile_args=extra_compile_args, extra_link_args=extra_link_args)



utils = ext(
            'utils.ext.libcontactutils',
            sources=['utils/ext/contactutils.cpp',
                     'utils/ext/boost_converters.cpp'
                     ]
        )

io = ext(
            'utils.ext.libio',
            sources=['utils/ext/io.cpp',
                     'utils/ext/boost_converters.cpp'
                     ]
        )

regularizer = ext(
            'coupling_prior.ext.libreg',
            sources=['coupling_prior/ext/Regularizer_PyWrapper.cpp',
                     'coupling_prior/ext/Regularizer.cpp',
                     'coupling_prior/ext/Parameters.cpp',
                     'utils/ext/boost_converters.cpp'
                     ]
        )

likelihood = ext(
            'coupling_prior.ext.libll',
            sources=['coupling_prior/ext/Likelihood_Dataset_PyWrapper.cpp',
                     'coupling_prior/ext/Likelihood_Dataset.cpp',
                     'coupling_prior/ext/Likelihood_Protein.cpp',
                     'coupling_prior/ext/Parameters.cpp',
                     'utils/ext/boost_converters.cpp'
                     ]
        )

counts = extc(
            'contact_prior.ext.counts.libmsacounts',
            sources=['contact_prior/ext/counts/msacounts.c']
        )


weighting = extc(
            'contact_prior.ext.weighting.libweighting',
            sources=['contact_prior/ext/weighting/weighting.c']
        )

setup(
    name="contact_prediction",
    version="1.0.0",
    description="cpp utils",
    license="AGPLv3",
    packages=find_packages(),
    ext_modules=[utils, io, regularizer, likelihood, counts, weighting],
    scripts=['contact_prior/grid_search_RF.py']
)
