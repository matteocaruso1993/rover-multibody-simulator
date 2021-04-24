# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 22:45:22 2021

@author: Matteo
"""

from setuptools import setup
from setuptools import find_packages




setup(name='rover_multibody_simulator',
      version='1.0.0',
      description='This package allows implementing a multibody dynamic simulator for a rover',
      url='https://github.com/matteocaruso1993/rover-multibody-simulator',
      author='Matteo Caruso',
      author_email='matteo.caruso@phd.units.it',
      license="GPLv3",
      packages=find_packages(),
      package_data = {'rover_multibody_simulator': ['data/config/config.ini','data/video/*.mp4','data/ground/*.json', 'data/model/*.pickle']},
      install_requires=['sympy==1.7.1','scipy==1.5.4', 'tqdm==4.51.0', 'numpy==1.19.2','matplotlib==3.3.3','dill==0.3.3', 'configparser==5.0.2',\
                        'cloudpickle==1.6.0','jsons==1.4.0','importlib==1.0.4','geneticalgorithm==1.0.2', 'numba==0.51.2', 'Cython==0.29.22'],
      zip_safe=False,
      include_package_data=True,
      python_requires='>=3.6'
      )
