from distutils.core import setup
setup(
  name = 'Pynoma',         
  packages = ['pynoma'],  
  version = '0.2.1',
  license='GNU General Public License v3.0',
  description = 'A Python API to communicate with gnomAD database.',   
  author = 'Felipe Colombelli, Paola Carneiro',                   
  author_email = 'bioinfo@hcpa.edu.br',
  url = 'https://github.com/bioinfo-hcpa/pynoma',
  keywords = ['gnomad', 'api', 'variants', 'genes'],
  install_requires=[
          'pandas>=1.0.5',
          'numpy>=1.19.0',
          'requests>=2.24.0',
          'seaborn'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3.8',
  ],
)

