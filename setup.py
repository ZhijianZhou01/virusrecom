# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r",encoding="utf-8") as fh:
  long_description = fh.read()

setuptools.setup(
  name="virusrecom",
  version="1.1.5",
  author="Zhi-Jian Zhou",
  author_email="zjzhou@hnu.edu.cn",
  description="An information-theory-based method for recombination detection of viral lineages.",
  keywords="recombination, virus, evolution, information entropy",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/ZhijianZhou01/virusrecom",
  packages=setuptools.find_packages(),
  install_requires=["matplotlib>=2.2.5",
                    "numpy>=1.19.3",
                    "pandas>=1.1.5",
                    "scipy>=1.5.4",
                    "psutil>=5.9.1",
                    "openpyxl>=3.0.5"
                     ],

  classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
        "Operating System :: OS Independent",
    ],
  entry_points={
             'console_scripts': [
                 'virusrecom = virusrecom.main:starts',
             ],
    }
)
