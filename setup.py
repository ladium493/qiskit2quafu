import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='qiskit2quafu',
    version='0.0.15',
    packages=setuptools.find_packages(),
    url='https://github.com/ladium493/qiskit2quafu',
    license='MIT',
    author='XR Wang',
    author_email='wang.x.aw@m.titech.ac.jp',
    description='This package converts the circuit of qiskit into a circuit supported by quafu',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['qiskit','pyquafu'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)