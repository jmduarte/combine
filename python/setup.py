from setuptools import setup, find_packages

setup(
    name="combine",
    version="1.0",
    description="Python analysis framework for CMS",
    long_description="Python package to facilitate analysis work in in the CMS collaboration.",
    url="http://github.com/guitargeek/geeksw",
    author="Jonas Rembser",
    author_email="jonas.rembser@cern.ch",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    # test_suite="tests",
    # setup_requires=["pytest-runner"],
    # tests_require=[],
    install_requires=[],
    scripts=['scripts/combineCards.py', 'scripts/commentUncerts.py', 'scripts/pruneUncerts.py', 'scripts/text2workspace.py', ],
)
