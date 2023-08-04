from setuptools import setup, find_packages
from setuptools import setup

# Reads the content of your README.md into a variable to be used in the setup below
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "rt", encoding="utf-8") as fh:
    install_requirements_txt = [line.strip() for line in fh.readlines()]

setup(
    name='sctools',
    version='0.0.0',    
    description='Single-cell utilities',
    url='https://github.com/ktroule/sctools',
    author='Kevin',
    author_email='nomail@nomail.com',
    license='',
    packages=['sctools', 'sctools.tools', 'sctools.metrics'],
    install_requires=install_requirements_txt,
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
