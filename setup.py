from os import path
from setuptools import setup,find_packages

exec(compile(open('gmgc_finder/gmgc_finder_version.py').read(),
             'gmgc_finder/gmgc_finder_version.py', 'exec'))


try:
    long_description = open('README.md', encoding='utf-8').read()
except:
    long_description = open('README.md').read()


# A bit hacky, but we want to have the output.md file in the package directory
# so it can be found by the setuptools machinery
if not path.exists('gmgc_finder/output.md'):
    from shutil import copyfile
    copyfile('docs/output.md', 'gmgc_finder/output.md')

setup(name='GMGC-Finder',
      version=__version__,
      description='Map genes and genome to the Global Microbial Gene Catalog (GMGC)',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
      ],
      url='https://github.com/psj1997/GMGC_Finder',
      author='Shaojun Pan',
      author_email='shaojun1997777@gmail.com',
      license='MIT',
      packages=['gmgc_finder'],
      install_requires=[
          'Biopython',
          'scikit-bio',
          'safeout',
          'tqdm',
      ],
      package_data={
             'gmgc_finder': ['*.md']},
      zip_safe=False,
      entry_points={
            'console_scripts': ['gmgc-finder=gmgc_finder.main:main'],
      }
      )
