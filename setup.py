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
      install_requires=['sympy','scipy', 'tqdm', 'numpy','matplotlib','dill', 'configparser', 'cloudpickle','jsons','importlib'],
      zip_safe=False,
      include_package_data=True,
      python_requires='>=3.6'
      )