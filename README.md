Combine
=======

## Introduction

This is a fork of the [HiggsAnalysis-CombinedLimit tool](https://github.com/cms-analysis/higgsanalysis-combinedlimit), developed by the CMS collaboration. The original tool is used by many analyses in the collaboration, not only in the Higgs sector. For more information how to use it, follow the documentation in the original repository.

![Brazilian plot](res/UpperLimit.png "Brazilian plot")

This fork was done on the 23. March 2019 from the `81x-root606` branch at the commit with the hash `5cc169efd9233011924e5fbd468ef05be044ed39`. The goal was:

* Easy standalone complilation with __cmake__
* modernize the tool by consistently applying the C++ 17 standard and porting all Python code to Python 3
* make the code nicer by applying __clang-format__ on the C++ code and formatting the Python code with __black__

With version 1.0 of this fork, the code differs from the original only by the points listed above plus the removal of non-essantial files from the repository. It will be decided in the future if compatibility with the original will be kept or if the functionality will diverge to ensure a modern workflow, flexibility and ease of use.


## Installation

Combine consists of both Python and C++ code, which have to be installed independently.

Two unusual build and runtime requirements, which are usually not installed in the usual linux system are:
* The [ROOT framework](https://root.cern.ch/)
* [CERN VDT](https://github.com/drbenmorgan/vdt) for fast math

For Arch Linux, these requirements are provided by the official repositories in the [root](https://www.archlinux.org/packages/community/x86_64/root/) and [cern-vdt](https://www.archlinux.org/packages/community/x86_64/cern-vdt/) packages.

### Python installation

You can install the python library and scripts just by executing the following command from the repository:
```
pip install --user python/
```

### C++ installation

The header files need to stay where they have been during the installation, so it's advised to copy the `libcombine` directory somewhere where you won't be bothered by it anymore:

```
cp -r libcombine $HOME/.local/combine
```

Now we can move to this directory and proceed with the build and installation:

```
cd $HOME/.local/combine
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local/ .
make -j8
make install
```

Lastly, make sure that your `LD_LIBRARY_PATH` contains `$HOME/.local/lib` and your `PATH` contains `$HOME/.local/bin` such that the binary and libraries can be found.

## Getting started

After installing the combine tool, you can test it by running the [Brazilian plot example](examples/brazilian_plots.py) to reproduce the plot at the top of this readme file. So far, this is the only supported example.

## For Developers

### Code formatting

To format the whole C++ code:
```
find libcombine -regex '.*\.\(cpp\|hpp\|cc\|cxx\|hh\|hxx\|h\)' -exec clang-format -style=file -i {} \;
```

To format the whole Python code:
```
black --line-length 120 python
```
