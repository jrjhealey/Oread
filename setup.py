from setuptools import setup
from oread.version import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='Oread',
    version=__version__,
    packages=['oread'],
    url='https://github.com/jrjhealey/Oread',
    license='GNU General Public License v3.0',
    author='Joe R. J. Healey',
    author_email='jrj.healey@gmail.com',
    long_description=long_description,
    long_description_content_type="text/markdown",
    description='A commandline and GUI tool for the creation of Artemis/Artemis Comparison Tool files.',
    install_requires=['biopython', 'gooey'],
    test_suite='nose.collector',
    tests_require=['nose'],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Environment :: MacOS X",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",

    ]

)
