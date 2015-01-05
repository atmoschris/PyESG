from setuptools import setup
import pyesg

setup(
    name = 'pyesg',
    version=pyesg.__version__,
    description = 'Python library of Earth Spherical Geometry (PyESG)',
    long_description=open('README.md').read(),
    author='Feng Zhu',
    author_email='feng.zhu@ssec.wisc.edu',
    url = 'https://github.com/lyricorpse/PyESG',
    license='BSD',
    py_modules = ['pyesg'],
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Programming Language :: Python :: 3'
    ]
)

