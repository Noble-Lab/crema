"""
Setup the crema package.
"""
import setuptools

with open("README.md", "r") as readme:
    LONG_DESC = readme.read()

DESC = (
    "Confidence estimation for peptide detection in mass spectrometry"
    " proteomics"
)

CATAGORIES = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setuptools.setup(
    name="crema-ms",
    author="William E. Fondrie",
    author_email="fondriew@gmail.com",
    description=DESC,
    long_description=LONG_DESC,
    long_description_content_type="text/markdown",
    url="https://github.com/Noble-Lab/crema",
    packages=setuptools.find_packages(),
    license="Apache 2.0",
    entry_points={"console_scripts": ["crema = crema.crema:main"]},
    classifiers=CATAGORIES,
    install_requires=["numpy>=1.18.1", "pandas>=1.0.3", "numba>=0.48.0"],
    use_scm_version=True,
    setup_requires=["setuptools-scm"],
    extras_require={
        "docs": [
            "numpydoc>=1.0.0",
            "sphinx-argparse>=0.2.5",
            "sphinx-rtd-theme>=0.5.0",
            "nbsphinx>=0.7.1",
            "ipykernel>=5.3.0",
        ],
        "dev": ["pre-commit>=2.7.1", "black>=19.10b0", "pytest>=6.0.1"],
    },
)
