from distutils.core import setup

setup(
    name='PyCyc',
    version='0.0.9',
    author='Eli Bogart',
    author_email='elb87@cornell.edu',
    packages=['pycyc','pycyc.test'],
    description='Python interface to Pathway Tools',
    long_description=open('README.txt').read()
    )
