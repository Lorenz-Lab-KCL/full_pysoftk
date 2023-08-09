from setuptools import setup, find_packages

setup(
    name='pySoftWhere',
    version='0.1.0',
    author='Robert M. Ziolek',
    author_email='robert.ziolek@kcl.ac.uk',
    packages= find_packages(exclude=["tests", "tests.*"]),
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    #url='http://pypi.python.org/pypi/TowelStuff/',
    #license='MIT',
    description='pySoftWhere: A Python Package for Analysing Interfaces of Soft Matter Nanostructures',
    #long_description=open('README.txt').read(),
    install_requires=[
        "networkx >= 2.6.3",
        "MDAnalysis >= 2.0.0",
        "scipy >= 1.7.2",
    ],
)
