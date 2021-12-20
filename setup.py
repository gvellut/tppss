from setuptools import find_packages, setup

with open("tppss/__init__.py") as f:
    for line in f:
        if line.find("__version__") >= 0:
            version = line.split("=")[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            break

with open("README.md") as f:
    readme = f.read()

with open("requirements.txt") as f:
    requirements = f.readlines()

with open("requirements-dev.txt") as f:
    requirements_dev = f.readlines()

with open("requirements-cli.txt") as f:
    requirements_cli = f.readlines()

setup_args = dict(
    name="tppss",
    version=version,
    description=("Compute sunrise / sunset times taking into account local topography"),
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/gvellut/tppss",
    author="Guilhem Vellut",
    author_email="g@vellut.com",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Utilities",
        "Topic :: Scientific/Engineering :: GIS",
    ],
    keywords="sunrise dem horizon",
    packages=find_packages(exclude=["docs", "tests"]),
    install_requires=requirements,
    extras_require={"dev": requirements_dev, "cli": requirements_cli},
    project_urls={
        "Bug Reports": "https://github.com/gvellut/tppss/issues",
        "Source": "https://github.com/gvellut/tppss",
    },
    entry_points={"console_scripts": ["tppss=tppss.cli:main"]},
)

setup(**setup_args)
