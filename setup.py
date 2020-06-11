from setuptools import setup,find_packages

exec(compile(open('gmgc_finder/gmgc_finder_version.py').read(),
             'gmgc_finder/gmgc_finder_version.py', 'exec'))

setup(name='GMGC-Finder',
      version=__version__,
      description='map genome to gmgc',
      long_description='',
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
             'docs': ['*.md']},
      zip_safe=False,
      entry_points={
            'console_scripts': ['gmgc-finder=gmgc_finder.main:main'],
      }
      )
