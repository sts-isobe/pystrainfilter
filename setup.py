from setuptools import setup

setup(
    name="pystrainfilter",
    version="0.1",
    entry_points={
        "console_scripts": [
            "pystrainfilter=pystrainfilter.pystrainfilter:main",
        ],
    },
    author="Michio Katouda",
    author_email="katouda@rist.or.jp",
    description="python interface of StrainFilter",
    url="https://github.com/mkatouda/pystrainfilter",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.8',
)
