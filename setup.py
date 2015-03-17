__author__ = 'Sander Bollen'

from setuptools import setup

setup(name="pyrefflat",
      version="0.1",
      description="A refFlat parser in pure python",
      author="Sander Bollen",
      author_email="sander@sndrtj.eu",
      license="MIT",
      packages=['pyrefflat'],
      zip_safe=False,
      scripts=['bin/bamCoverage', 'bin/gVCFCoverage', 'bin/createMargin', 'bin/refFlat2Bed'])