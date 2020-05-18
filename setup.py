from setuptools import setup

from setuptools import setup

setup(name='Genome2gmgc',
      version='0.1',
      description='map genome to gmgc',
      long_description='',
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
      ],
      url='https://github.com/psj1997/Genome2gmgc',
      author='Shaojun Pan',
      author_email='shaojun1997777@gmail.com',
      license='MIT',
      packages=['genome2gmgc'],
      install_requires=[
          'Biopython',
          'scikit-bio',
            'tqdm'
      ],
      include_package_data=True,
      zip_safe=False,
      entry_points={
            'console_scripts': ['genome2gmgc=genome2gmgc.main:main'],
      }
      )