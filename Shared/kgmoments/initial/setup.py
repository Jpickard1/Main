from setuptools import setup, Extension

setup(name = "Krongraph",
      version = "1.0",
      ext_modules = [Extension("krongraph", ["krongraph.cc"])])
      
      
