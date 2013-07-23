#!/usr/bin/env python

from distutils.core import setup

setup(name='networkb',
      version='0.1dev0',
      description='Network analysis for brain data',
      author='Andres Babino',
      author_email='ababino@gmail.com',
      packages=['networkb', 'networkb.algorithms','networkb.classes','networkb.ploters'],
     )
