from setuptools import setup, find_packages

setup(
      name		='tissuespecific',
      version		='0.3',
      description	='Reconstruction and analysis of tissue-specific metabolic models',
      classifiers	=['Development Status :: 3 - Alpha',
        		  'License :: OSI Approved :: MIT License',
        		  'Programming Language :: Python :: 3.6',
        		  'Topic :: Metabolic-modelling :: Metabolic-engineering'],
      url		='https://github.com/Andr3aCbb/tissuespecific',
      author		='Andrea Cabbia',
      author_email	='a.cabbia@tue.nl',
      license		='MIT',
      packages		= find_packages(),
      install_requires	=['cobra', 'pandas', 'GEOparse', 'cobrababel'],
      zip_safe		=False
)



