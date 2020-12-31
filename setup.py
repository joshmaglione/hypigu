from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="hypigu",
    version="1.0",
    description="A SageMath package that provides functions to compute the Igusa local zeta function associated with hyperplane arrangements.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/joshmaglione/hypigu",
    author="Joshua Maglione",
    author_email="joshmaglione@gmail.com",
    license="MIT",
    classifiers=[
        'Topic :: Mathematics',
        "License :: MIT License",
        "Programming Language :: Python 3",
    ],
    packages=["hypigu"],
    include_package_data=True,
)
