from setuptools import setup, find_packages

setup(
    name='sdi_variation',
    version='0.1.0',
    packages=find_packages(),  # automatically finds `your_module` and submodules
    install_requires=[],       # add any dependencies here
    author='Isaac von Riedemann',
    author_email='imv@sfu.ca',
    description='Code to regenerate analysis and utility functions for computing annotator variability on spatialLIBD',
    long_description=open('../README.md').read(),  # optional
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)