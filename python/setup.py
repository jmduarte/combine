from setuptools import setup, find_packages

setup(
    name="combine",
    version="1.0",
    description="CMS Higgs Combination toolkit.",
    long_description="CMS Higgs Combination toolkit.",
    url="http://github.com/guitargeek/combine",
    author="The CMS collaboration",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=[],
    scripts=['scripts/combineCards.py', 'scripts/commentUncerts.py', 'scripts/pruneUncerts.py', 'scripts/text2workspace.py', ],
)
