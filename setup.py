from setuptools import setup, find_packages

setup(
    name='logisticsusie',
    version='1.0.0',
    url='https://github.com/karltayeb/logisticsusie.git',
    author='Karl Tayeb',
    author_email='karl.tayeb7@gmail.com',
    description='Python implimentation of logistic susie',
    packages=find_packages(),    
    install_requires=['numpy >= 1.11.1', 'matplotlib >= 1.5.1'],
)