from distutils.core import setup
setup(
  name = 'pynomad',         
  packages = ['pynomad'],  
  version = '0.0.1',
  license='GNU General Public License v3.0',
  description = 'A Python API to communicate with gnomAD database.',   
  author = 'Felipe Colombelli',                   
  author_email = 'fcolombelli@inf.ufrgs.br',
  url = 'https://github.com/bioinfo-hcpa/pynomad',
  keywords = ['gnomad', 'api', 'variants', 'genes'],
  install_requires=[
          'pandas>=1.0.5',
          'numpy>=1.19.0',
          'requests>=2.24.0'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: GNU General Public License v3.0',
    'Programming Language :: Python :: 3.8',
  ],
)

