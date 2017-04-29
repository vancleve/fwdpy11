from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os,glob

__version__ = '0.1.1a0'

#clang/llvm is default for OS X builds.
#can over-ride darwin-specific options
#with setup.py --gcc install
if '--gcc' in sys.argv:
    USE_GCC = True
    sys.argv.remove('--gcc')
else:
    USE_GCC = False

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

PKGS=['fwdpy11']

INCLUDES=[
    'fwdpy11/headers',
    'fwdpy11/headers/fwdpp',
    # Path to pybind11 headers
    get_pybind_include(),
    get_pybind_include(user=True),
    os.path.join(sys.prefix, 'include')
]

LIBRARY_DIRS=[
    os.path.join(sys.prefix, 'lib')
    ]
ext_modules = [
    Extension(
        'fwdpy11.fwdpp_types',
        ['fwdpy11/src/fwdpp_types.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.fwdpp_extensions',
        ['fwdpy11/src/fwdpp_extensions.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.fwdpy11_types',
        ['fwdpy11/src/fwdpy11_types.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.sampling',
        ['fwdpy11/src/fwdpy11_sampling.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.fitness',
        ['fwdpy11/src/fwdpy11_fitness.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.trait_values',
        ['fwdpy11/src/fwdpy11_trait_values.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
        Extension(
        'fwdpy11.wfevolve',
        ['fwdpy11/src/wfevolve.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.wfevolve_qtrait',
        ['fwdpy11/src/wfevolve_qtrait.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.gsl_random',
        ['fwdpy11/src/gsl_random.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    Extension(
        'fwdpy11.multilocus',
        ['fwdpy11/src/multilocus.cc'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    ),
    ]
    


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin' and USE_GCC is False:
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden') and (sys.platform != 'darwin' or USE_GCC is True):
                opts.append('-fvisibility=hidden')
            if has_flag(self.compiler,'-g0'):
                opts.append('-g0')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
            if sys.platform == 'darwin' and USE_GCC is False:
                ext.extra_link_args=['-stdlib=libc++', '-mmacosx-version-min=10.7']
        build_ext.build_extensions(self)

#Figure out the headers we need to install:
generated_package_data={}
for root, dirnames, filenames in os.walk('fwdpy11/headers'):
    if 'testsuite' not in root and 'examples' not in root and 'python_examples' not in root:
        g=glob.glob(root+'/*.hh')
        if len(g)>0:
            replace=root.replace('/','.')
            if replace not in PKGS:
                PKGS.append(replace) #If there's a header file, we add the directory as a package
            generated_package_data[replace]=['*.hh']
        g=glob.glob(root+'/*.hpp')
        if len(g)>0:
            replace=root.replace('/','.')
            if replace not in PKGS:
                PKGS.append(replace) #If there's a header file, we add the directory as a package
            try:
                if '*.hpp' not in generated_package_data[replace]:
                    generated_package_data[replace].append('*.hpp')
            except:
                generated_package_data[replace]=['*.hpp']
        g=glob.glob(root+'/*.tcc')
        if len(g)>0:
            replace=root.replace('/','.')
            if replace not in PKGS:
                PKGS.append(replace) #If there's a template implementation file, we add the directory as a package
            try:
                if '*.tcc' not in generated_package_data[replace]:
                    generated_package_data[replace].append('*.tcc')
            except:
                generated_package_data[replace]=['*.tcc']

long_desc = open("README.rst").read()

setup(
    name='fwdpy11',
    version=__version__,
    author='Kevin Thornton',
    author_email='krthornt@uci.edu',
    url='http://molpopgen.github.io/fwdpy11',
    classifiers=['Intended Audience :: Science/Research',
               'Topic :: Scientific/Engineering :: Bio-Informatics',
               'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'],
    description='Forward-time population genetic simulation in Python',
    license='GPL >= 3',
    requires=['pybind11','numpy'],
    provides=['fwdpy11'],
    obsoletes=['none'],
    data_files=[('fwdpy11',['COPYING', 'README.rst'])],
    long_description=long_desc,
    ext_modules=ext_modules,
    install_requires=['pybind11>=1.7'],
    cmdclass={'build_ext': BuildExt},
    packages=PKGS,
    package_data=generated_package_data,
    zip_safe=False,
)
