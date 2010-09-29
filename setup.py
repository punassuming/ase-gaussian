from setuptools import setup, find_packages
import sys, os

version = '0.4'

setup(name='rca',
      version=version,
      description="Computational Wrapper for Gaussian and Jacapo Calculations",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='Rich Alesi',
      author_email='walesi at cmu.edu',
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
