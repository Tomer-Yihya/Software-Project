from setuptools import Extension, setup

module = Extension("mySpkmeansModule",sources =['spkmeans.c','spkmeansmodule.c'])
setup(name='mySpkmeansModule',
      version='1.0',
      description='Python wrapper for C extension',
      ext_modules=[module])

