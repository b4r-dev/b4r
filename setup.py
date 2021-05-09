from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# long_description(後述)に、GitHub用のREADME.mdを指定
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='b4rpipe', # パッケージ名(プロジェクト名)
    packages=['b4rpipe'], # パッケージ内(プロジェクト内)のパッケージ名をリスト形式で指定

    version='0.0.10', # バージョン

    license='MIT', # ライセンス

    install_requires=['numpy','scipy','astropy','xarray','matplotlib','fmflow','bokeh','netcdf4','python-casacore','decode','tqdm'], # pip installする際に同時にインストールされるパッケージ名をリスト形式で指定

    author='astroysmr', # パッケージ作者の名前
    author_email='astro.yoshimura@gmail.com', # パッケージ作者の連絡先メールアドレス

    url='https://github.com/LMT-heterodyne-dev/b4rpipe', # パッケージに関連するサイトのURL(GitHubなど)

    description='Pipeline reduction tools for B4R (2mm heterodyne reciver) on LMT 50m @Mexico', # パッケージの簡単な説明
    long_description=long_description, # PyPIに'Project description'として表示されるパッケージの説明文
    long_description_content_type='text/markdown', # long_descriptionの形式を'text/plain', 'text/x-rst', 'text/markdown'のいずれかから指定
    keywords='B4R LMT', # PyPIでの検索用キーワードをスペース区切りで指定

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
    ],
)
