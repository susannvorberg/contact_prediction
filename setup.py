
import numpy.distutils.intelccompiler
from setuptools import setup, Extension, find_packages


def ext(name, sources=[],
        include_dirs=[],
        library_dirs=[],
        libraries=[],
        extra_compile_args=['-Wall -g -shared -O2 -fPIC -fopenmp -Wl,--export-dynamic -std=c++11'],
        extra_link_args=['-lpython2.7 -lboost_python -lboost_iostreams  -lpthread -larmadillo -lmsgpack -fopenmp']):
    return Extension(name,
                     include_dirs=include_dirs,
                     library_dirs=library_dirs,
                     libraries=libraries,
                     sources=sources,
                     extra_compile_args=extra_compile_args,
                     extra_link_args=extra_link_args)


setup(
    name="Contact Prediction Extension Modules",
    version="1.0.0",
    description="bayes_utils_py",
    license="AGPLv3",
    packages=find_packages(),
    ext_modules=[
        ext(
            'build.libcontactutils',
            sources=['cpp_modules/contactutils.cpp',
                     'cpp_modules/util_math.cpp',
                     'cpp_modules/boost_converters.cpp']
        )
    ],
    scripts=['benchmark/append_to_evaluation_file.py']
)
