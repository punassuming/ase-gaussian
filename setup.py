from setuptools import setup, find_packages
import sys, os

version = '0.15'

setup(name='gauss',
      version=version,
      description="gaussian parser for integration with ASE objects",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='qchem gaussian parse ase',
      author='Rich Alesi',
      author_email='walesi',
      url='',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
