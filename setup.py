import re
import sys

try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit(
        'We need the Python library setuptools to be installed. '
        'Try runnning: python -m ensurepip'
    )

with open("README.md", 'r') as fin:
    long_description = fin.read()

with open("README.md", 'r') as fin:
    version_line_regex = re.compile(r'^\s*__version__\s*=\s*[\'"]([^\'"]+)[\'"]')
    for line in fin:
        match = version_line_regex.match(line)
        if match:
            version = match.group(1)

setup(
    name="jasentool",
    version=version,
    description="A mongodb validation tool for comaparing pipeline outputs",
    long_description_markdown_filename=long_description,
    long_description_content_type="text/markdown",
    author="Ryan James Kennedy",
    license="GPLv3",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.10",
    ],
    install_requires=["pymongo", "openpyxl", "biopython"],
    entry_points={"console_scripts": ["jasentool=jasentool.__main__:main"]},
    packages=find_packages(exclude=("tests")),
)