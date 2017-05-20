
import numpy.distutils.intelccompiler
from setuptools import setup, Extension, find_packages



#-D_GLIBCXX_USE_CXX11_ABI=0: because GCC 5 issue with dual ABI: https://gcc.gnu.org/onlinedocs/gcc-5.2.0/libstdc++/manual/manual/using_dual_abi.html

def ext(name,
        sources=[],
        include_dirs=['cpp_modules', '/home/vorberg/anaconda2/envs/py27/include/python2.7', '/home/vorberg/anaconda2/envs/py27/include', '/usr/local/include'],
        library_dirs=['/home/vorberg/anaconda2/envs/py27/lib'],
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

module1 = ext(
            'build.libcontactutils',
            sources=['cpp_modules/contactutils.cpp',
                     'cpp_modules/boost_converters.cpp'
                     ]
        )

module2 = ext(
            'build.libio',
            sources=['cpp_modules/io.cpp',
                     'cpp_modules/boost_converters.cpp'
                     ]
        )

module3 = ext(
            'coupling_prior.ext.libreg',
            sources=['coupling_prior/ext/Regularizer_PyWrapper.cpp',
                     'coupling_prior/ext/Regularizer.cpp',
                     'coupling_prior/ext/Parameters.cpp',
                     'cpp_modules/boost_converters.cpp'
                     ]
        )

module4 = ext(
            'coupling_prior.ext.libll',
            sources=['coupling_prior/ext/Likelihood_Dataset_PyWrapper.cpp',
                     'coupling_prior/ext/Likelihood_Dataset.cpp',
                     'coupling_prior/ext/Likelihood_Protein.cpp',
                     'coupling_prior/ext/Parameters.cpp',
                     'cpp_modules/boost_converters.cpp'
                     ]
        )




setup(
    name="Contact Prediction Extension Modules",
    version="1.0.0",
    description="cpp utils",
    license="AGPLv3",
    packages=find_packages(),
    ext_modules=[module1, module2, module3, module4]
)
