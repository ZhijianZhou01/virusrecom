import setuptools

with open("README.md", "r",encoding="utf-8") as fh:
  long_description = fh.read()

setuptools.setup(
  name="virusrecom",
  version="1.1.0.0",
  author="Zhi-Jian Zhou",
  author_email="zjzhou@hnu.edu.cn",
  description="An information-theory-based method for recombination detection of viral lineages.",
  keywords='recombination', 
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/ZhijianZhou01/virusrecom",
  packages=setuptools.find_packages(),
  install_requires=["matplotlib","pandas","numpy", "scipy"],
  include_package_data=True,
  
  classifiers=[
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
