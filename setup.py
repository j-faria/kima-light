#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='pykimalight',
      version=open('VERSION').read().strip(), # same as kima
      description='Analysis of results from kima',
      author='João Faria',
      author_email='joao.faria@astro.up.pt',
      license='MIT',
      url='https://github.com/j-faria/kima-light/tree/master/pykima',
      packages=['pykimalight'],
      install_requires=[
        'numpy',
        'scipy',
        'matplotlib>=1.5.3',
        'corner',
      ],
      entry_points={
        'console_scripts': [
            'kimal-showresults = pykimalight.showresults:showresults',
            'kimal-checkpriors = pykimalight.check_priors:main',
            'kimal-template = pykimalight.make_template:main',
            ]
        },
     )
