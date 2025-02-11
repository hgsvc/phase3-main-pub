from setuptools import setup, find_packages

setup(
    name='MHC-Annotation',
    version='0.1',
    license='MIT',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Torsten Houwaart',
    author_email='houwaart@hhu.de',
    url='https://github.com/DiltheyLab/MHC-annotation',
    packages=find_packages(),
    package_dir={'mhca': 'mhca'},
    package_data={'': ['data/*']},
    entry_points={
        'console_scripts': [
            'mhca=mhca.__main__:main',
        ],
    },
    python_requires=">=3.7",
)
