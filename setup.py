import os
import subprocess
import sys
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

module_name = 'cosmicRays'
# Command line flags forwarded to CMake (for debug purpose)
cmake_cmd_args = []
for f in sys.argv:
    if f.startswith('-D'):
        cmake_cmd_args.append(f)

for f in cmake_cmd_args:
    sys.argv.remove(f)

def _get_env_variable(name, default='OFF'):
    if name not in os.environ.keys():
        return default
    return os.environ[name]

class CMakeExtension(Extension):
    def __init__(self, name, cmake_lists_dir='.', **kwarg):
        Extension.__init__(self, name, sources=[], **kwarg)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)

class CMakeBuild(build_ext):
    def build_extensions(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("Cannot find cmake")

        for ext in self.extensions:
            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cfg = 'Debug' if _get_env_variable('COSMICRAYS_DEBUG') == 'ON' else 'Release'
            cmake_args = [
                '-DCMAKE_BUILD_TYPE=%s' % cfg,
                # Ask CMake to place the resulting library in the directory
                # containing the extension
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                # Other intermediate static libraries are placed in a
                # temporary build directory instead
                '-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), self.build_temp),
                # Hint CMake to use the same Python executable that
                # is launching the build, prevents possible mismatching if
                # multiple versions of Python are installed
                '-DPYTHON_EXECUTABLE={}'.format(sys.executable),
                # Add other project-specific CMake arguments if needed
                # ...
            ]

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            # Config
            subprocess.check_call(['cmake', ext.cmake_lists_dir] + cmake_args,
                                  cwd=self.build_temp)
            # Build
            subprocess.check_call(['cmake', '--build', '.', '--config', cfg],
                                  cwd=self.build_temp)

version = "0.1.0"

setup(name='cosmicRays',
      packages=['cosmicRays'],
      version=version,
      description='Port of the LSST cosmic ray finder code',
      author='Douglas N Friedel',
      install_requires=[],
      ext_modules=[CMakeExtension(module_name)],
      cmdclass={'build_ext': CMakeBuild},
      package_dir={"", "python"},
      )
