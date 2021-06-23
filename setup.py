import setuptools
from io import open

requirements = [
    'numpy',
    'pytest',
    'setuptools'
]

setuptools.setup(name='githubActions',
      maintainer='Chunhong Yoon',
      version='0.0.1',
      maintainer_email='yoon82@stanford.edu',
      description='Python application',
      long_description=open('README.md', encoding='utf8').read(),
      url='https://github.com/chuckie82/githubActions.git',
      packages=setuptools.find_packages(),
      install_requires=requirements,
      classifiers=[
          "Programming Language :: Python :: 3",
          "Operating System :: OS Independent",
      ],
      zip_safe=False)
