import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="geneannotator-burgoa",
    version="0.0.1",
    author="Javier Burgoa",
    author_email="burgoacardasjavier@gmail.com",
    description="A package for functional annotation of genes based on available public resources like emapper and EggNOG",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.embl.de/burgoa/geneannotator.git",
    project_urls={
        "Bug Tracker": "https://git.embl.de/burgoa/geneannotator/-/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
