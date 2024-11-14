"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='arteryfe',  # Required
    version='1.0',  # Required
    description='Implementation of the 1D blood flow equations in FEniCSx.',  # Optional
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',  # Optional (see note above)
    url='https://github.com/yourusername/arteryfe',  # Optional
    author='Your Name',  # Optional
    author_email='your.email@example.com',  # Optional
    classifiers=[  # Optional
        'Development Status: 4 - Beta',
        'Intended Audience: Developers',
        'Topic: Scientific/Engineering',
        'License: OSI Approved :: MIT License',
        'Programming Language: Python :: 3',
        'Programming Language: Python :: 3.7',
        'Programming Language: Python :: 3.8',
        'Programming Language: Python :: 3.9',
    ],
    keywords='fenicsx blood flow simulation',  # Optional
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required
    python_requires='>=3.7, <4',
    install_requires=['numpy', 'mpi4py','petsc4py' , 'matplotlib', 'scipy'],  # Optional
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/yourusername/arteryfe/issues',
        'Source': 'https://github.com/yourusername/arteryfe',
    },
)