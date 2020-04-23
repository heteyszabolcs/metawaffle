import os
from setuptools import setup, find_packages, Command


__version__ = None
exec(open('metawaffle/version.py').read())

class CleanCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info ./htmlcov')

setup(
    name='metawaffle',
    version=__version__,
    description='Deconvolution of regions of interest according to their structurual pattern.',
    url='https://github.com/3DGenomes/meta-waffle',
    setup_requires=[
        'setuptools>=18.0',
    ],
    packages=find_packages(),
    install_requires=[
        'cython',
        'numpy==1.16.3',
        'scikit-learn',
        'neupy',
        'matplotlib',
        'pysam',
        'scipy',
        'pandas',
        'pickle-mixin',
        'functools',
        'multiprocess',
    ],
    scripts=['bin/metawaffle'],
    cmdclass={
        'clean': CleanCommand
    },
    classifiers=(
        "Programming Language :: Python :: 2.7",
    ),
)
