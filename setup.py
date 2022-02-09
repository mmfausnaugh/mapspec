
from setuptools import setup

setup(
    name='mapspec',
    version='0.2.0',
    author='Michael Fausnaugh',
    author_email='faus@mit.edu',
    url='https://github.com/mmfausnaugh/mapspec',
    license='MIT',
    packages=['mapspec'],
    scripts=['scripts/make_ref',
             'scripts/smooth_ref',
             'scripts/do_map'],
#    package_data = ['examples/test_data/', 'examples/mapspec_test' ]

    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'matplotlib',
        
    ]
)
