from os import path
from setuptools import setup,find_packages

exec(compile(open('gmgc_mapper/gmgc_mapper_version.py').read(),
             'gmgc_mapper/gmgc_mapper_version.py', 'exec'))


try:
    long_description = open('README.md', encoding='utf-8').read()
except:
    long_description = open('README.md').read()


setup(name='GMGC-mapper',
      version=__version__,
      description='Map genes and genome to the Global Microbial Gene Catalog (GMGC)',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
      ],
      url='https://github.com/BigDataBiology/GMGC-mapper',
      author='Shaojun Pan',
      author_email='shaojun1997777@gmail.com',
      license='MIT',
      packages=['gmgc_mapper'],
      install_requires=[
          # Technically, numpy is not directly needed, but some downstream
          # dependencies use it and fail to declare they need it:
          'numpy',
          'Biopython',
          'scikit-bio',
          'tqdm',
          'pyyaml',
          'atomicwrites',
      ],
      package_data={
             'gmgc_mapper': ['*.md']},
      zip_safe=False,
      entry_points={
            'console_scripts': ['gmgc-mapper=gmgc_mapper.main:main'],
      }
      )
