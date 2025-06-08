from setuptools import setup, find_packages

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="nudca",
    version="0.1.0",
    description="A Numerical Code for Nuclear Decay Chains in Astrophysics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Qiuhong-Chen",
    author_email="chohonche@163.com",
    url="https://github.com/nudca-code/NuDCA",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "numba"
    ],
    include_package_data=True,
    package_data={
        "nudca": [
            "data/*.npz", "data/*.json", "data/*.csv", "data/*.xlsx"
        ],

        "nudca.kilonovae": [
            "Tanaka_Opacity_Data/*.txt",
            "Tanaka_Opacity_Data/readme.txt",
            "Tanaka_Opacity_Data/Ye0.10/*.txt",
            "Tanaka_Opacity_Data/Ye0.15/*.txt",
            "Tanaka_Opacity_Data/Ye0.20/*.txt",
            "Tanaka_Opacity_Data/Ye0.25/*.txt",
            "Tanaka_Opacity_Data/Ye0.30/*.txt",
            "Tanaka_Opacity_Data/Ye0.35/*.txt",
            "Tanaka_Opacity_Data/Ye0.40/*.txt",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    keywords="nuclear decay astrophysics kilonova r-process",
    project_urls={
        "Documentation": "https://nudca-code.github.io/nudca.github.io/",
        "Source": "https://github.com/nudca-code/NuDCA",
        "Tracker": "https://github.com/nudca-code/NuDCA/issues",
    },
)
