#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='networkb',
      version='0.1dev0',
      description='Network analysis for brain data',
      author='Andres Babino',
      author_email='ababino@gmail.com',
      packages=find_packages(),
      include_package_data=True,
      exclude_package_data = {'': ['.gitignore']},
     )
