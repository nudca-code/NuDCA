#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

# Read README.md for long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements from pyproject.toml
def read_requirements():
    with open("pyproject.toml", "r", encoding="utf-8") as f:
        content = f.read()
        # Extract dependencies from pyproject.toml
        start = content.find("[tool.poetry.dependencies]")
        end = content.find("[tool.poetry.group.dev.dependencies]")
        if start == -1 or end == -1:
            return []
        
        deps = content[start:end].split("\n")
        requirements = []
        for dep in deps:
            if "=" in dep and not dep.startswith("#"):
                package = dep.split("=")[0].strip()
                version = dep.split("=")[1].strip().strip('"').strip("'")
                if package != "python":
                    requirements.append(f"{package}{version}")
        return requirements

setup(
    name="nudca",
    version="0.1.0",
    author="Qiuhong-Chen",
    author_email="chohonche@163.com",
    description="A Numerical Code for Nuclear Decay Chains in Astrophysics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chohonche/nudca",
    project_urls={
        "Documentation": "https://nudca.readthedocs.io/",
        "Source": "https://github.com/chohonche/nudca",
        "Tracker": "https://github.com/chohonche/nudca/issues",
    },
    packages=find_packages(include=["nudca", "nudca.*"]),
    package_data={
        "nudca": [
            "data/*.csv",
            "data/*.json",
            "data/*.h5",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Nuclear Physics",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.2.4",
            "pytest-cov>=2.12.0",
            "flake8>=3.9.2",
            "black>=22.3.0",
            "isort>=5.10.0",
            "mypy>=0.910",
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "nudca=nudca.cli:main",
        ],
    },
    zip_safe=False,
)
