#from distutils.core import setup
from setuptools import setup

setup(
    name = 'pyesg',
    py_modules = ['pyesg'],
    version = '0.0.1',
    description = 'Python library of Earth Spherical Geometry (PyESG)',
    author='Feng Zhu',
    author_email='feng.zhu@ssec.wisc.edu',
    url = 'https://github.com/lyricorpse/PyESG',
    license='ISC',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Programming Language :: Python :: 3'
    ]
)

