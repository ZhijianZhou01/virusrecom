import setuptools

with open("README.md", "r",encoding="utf-8") as fh:
  long_description = fh.read()

setuptools.setup(
  name="virusrecom",
  packages=["virusrecom"],
  install_requires=["matplotlib","pandas","numpy", "scipy", "argparse", "psutil","platform"],
  version="1.0",
  author="Zhi-Jian Zhou",
  author_email="zjzhou@hnu.edu.cn",
  description="An information-theory-based method for recombination detection of viral lineages.",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/ZhijianZhou01/virusrecom",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3.7",
  "License :: OSI Approved :: GNU General Public License v3.0",
  "Operating System :: OS Independent",
  ],
)