import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    install_requires = [line.strip() for line in fh.readlines() if line.strip()]

exec(open('atomict/version.py').read())

setuptools.setup(
    name="atomict",
    version=__version__,
    author="Alain Richardt",
    author_email="alain@atomictessellator.com",
    description="The client application for the https://atomictessellator.com/",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/your-cli-project",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'at=atomict.cli:main',
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires='>=3.6',
    install_requires=install_requires,
)