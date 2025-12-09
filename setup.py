import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="aKmerBroom",
    version="1.0.0",
    author="Camila Duitama González",
    author_email="cduitama@pasteur.fr",
    description="Ancient oral DNA decontamination using Bloom filters on k-mer sets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CamilaDuitama/aKmerBroom",
    project_urls={
        "Bug Tracker": "https://github.com/CamilaDuitama/aKmerBroom/issues",
        "Paper": "https://www.cell.com/iscience/pdf/S2589-0042(23)02134-X.pdf"
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "biopython>=1.78",
        "pybloomfiltermmap3>=0.5.0",
    ],
    entry_points={
        "console_scripts": [
            "aKmerBroom=akmerbroom:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt", "*.yml"],
    },
)