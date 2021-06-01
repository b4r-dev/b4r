from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# long_description(後述)に、GitHub用のREADME.mdを指定
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='b4rpipe',
    packages=['b4rpipe'],

    version='0.1.2', 

    license='MIT',

    install_requires=['numpy','scipy','astropy','xarray','matplotlib','fmflow','bokeh','netcdf4','python-casacore','decode','tqdm'], # pip installする際に同時にインストールされるパッケージ名をリスト形式で指定

    author='astroysmr',
    author_email='astro.yoshimura@gmail.com',

    url='https://github.com/LMT-heterodyne-dev/b4rpipe',

    description='Pipeline reduction tools for B4R (2mm heterodyne reciver) on LMT 50m @Mexico',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='B4R LMT',

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],
)
