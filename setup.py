import setuptools
import platform
import subprocess

from setuptools.command.build import build
from setuptools import  find_namespace_packages


'''
Template from https://python.plainenglish.io/building-hybrid-python-c-packages-8985fa1c5b1d
based an modification of  https://www.benjack.io/2017/06/12/python-cpp-tests.html
'''

libname='qoptkit'

class CustomBuild(build):
    def run(self):
        try:
            if platform.system() == "Darwin":
                subprocess.check_output(['clang', '--version'])
            else:
                subprocess.check_output(['g++', '--version'])
        except OSError:
            raise RuntimeError("g++ must be installed to build qoptkit." )
            
        try:
            subprocess.check_output(['ar', '--version'])
        except OSError:
            raise RuntimeError("ar tool must be installed to build qoptkit." )
        
        try:
            subprocess.check_output(['make', '--version'])
        except OSError:
            raise RuntimeError("Make must be installed to build qoptkit." )
        

        subprocess.check_call(['make'],cwd='./src/' + libname )


setuptools.setup(
       cmdclass={
            "build": CustomBuild,
        },
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        libname : ["*"],
    }
)