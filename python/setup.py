from codecs import open as codecs_open
from setuptools import setup, find_packages




setup(name='BigVAR',
      version='0.0.1',
      description=u"BigVAR Python Port",
      install_requires=['numpy','statsmodels','numba'],
      classifiers=[],
      keywords='',
      author=u"Will Nicholson",
      author_email='wbn8@cornell.edu',
      url='https://github.com/wbnicholson/BigVAR',
      license='GPL>=2',
      packages=['BigVAR'],
      
      )
