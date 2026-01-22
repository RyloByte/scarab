import re
from pathlib import Path

from setuptools import setup

with open("README.md", "r") as readme:
    LONG_DESCRIPTION = readme.read()

CLASSIFIERS = [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

pks = ['scarab']


def read_version():
    version_file = Path(__file__).parent / "src" / "scarab" / "__init__.py"
    contents = version_file.read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*[\'"]([^\'"]+)[\'"]', contents, re.M)
    if not match:
        raise RuntimeError("Unable to find version in src/scarab/__init__.py")
    return match.group(1)


SETUP_METADATA = \
    {
        "name": "scarab",
        "version": read_version(),
        "description": "Software for recruiting metagenomic reads using single-cell amplified genome data.",
        "long_description": LONG_DESCRIPTION,
        "long_description_content_type": "text/markdown",
        "author": "Ryan McLaughlin, Connor Morgan-Lang",
        "author_email": "mclaughlinr2@gmail.com",
        "url": "https://github.com/hallamlab/scarab",
        "license": "GPL-3.0",
        "python_requires": ">=3.8",
        "include_package_data": True,
        "package_dir": {'': 'src'},  # Necessary for proper importing
        "packages": pks,
        "package_data": {
            'scarab': ['configs/*']},
        "entry_points": {'console_scripts': ['scarab = scarab.__main__:main']},
        "classifiers": CLASSIFIERS,
        "install_requires": [
            "numpy",
            "pandas",
            "scikit-learn",
            "umap-learn",
            "hdbscan",
            "sourmash",
            "pyfastx",
            "scikit-bio",
            "tqdm",
        ]
    }

setup(**SETUP_METADATA)
