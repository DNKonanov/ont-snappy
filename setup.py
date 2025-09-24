import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ont-snappy",
    version="0.1.5",
    author="D.N. Konanov",
    author_email="konanovdmitriy@gmail.com",
    description="Nanopore-based methylation sites caller",
    long_description="snappy",
    long_description_content_type="",
    url="https://github.com/DNKonanov/ont-snappy",
    project_urls={
        "Bug Tracker": "https://github.com/DNKonanov/ont-snappy",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    include_package_data=True,
    packages=['snappy', 'snappy.src'],
    install_requires=[
        'biopython',
        'matplotlib',
        'scipy',
        'seaborn',
        'tqdm',
        'polars',
        'pandas'

    ],
    entry_points={
        'console_scripts': [
            'snappy=snappy.snappy:main'
        ]
    }
)
