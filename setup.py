import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

nocpp = False

includepath = os.getenv("COMBINE_INCLUDE_DIR")

if os.getenv("COMBINE_LIB_DIR") is None and includepath is None:
    raise ValueError("Environment variable COMBINE_INCLUDE_PATH not set! Please set it to the permanent location of the combine header files.")

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DCMAKE_INSTALL_LIBDIR=" + extdir,
            "-DCOMBINE_HEADERS_DIR=" + includepath,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
        build_args += ["--", "-j8"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get("CXXFLAGS", ""), self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

ext_modules = [CMakeExtension("_combine")]

if not os.getenv("COMBINE_LIB_DIR") is None:
    ext_modules.pop(0)

setup(
    name="combine",
    version="1.0",
    description="CMS Higgs Combination toolkit.",
    long_description="CMS Higgs Combination toolkit.",
    url="http://github.com/guitargeek/combine",
    author="The CMS collaboration",
    packages=find_packages(),
    include_package_data=True,
    ext_modules=ext_modules,
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=[],
    scripts=[
        "scripts/combineCards.py",
        "scripts/commentUncerts.py",
        "scripts/pruneUncerts.py",
        "scripts/text2workspace.py",
    ],
)
