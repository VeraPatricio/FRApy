import setuptools

setuptools.setup(
    name="FRApy",
    version="0.0.3",
    author="Vera Patricio",
    author_email="verapatricio90@gmail.com",
    description="FRApy: Fitting Resolved Arcs with python",
    long_description="FRApy is python code designed to fit gravitational arcs in image plane with analytical modes, such as metallicity gradients and velocity fields, taking into account the lensing distortions.",
    long_description_content_type="text/markdown",
    url=[
        "https://github.com/VeraPatricio/FRApy",
        "https://readthedocs.org/projects/frapy/"
         ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['numpy',
                        'astropy',
                        'reproject',
                        'emcee',
                        'pickle',
                        'matplotlib',
                        'corner',
                        ],
)
