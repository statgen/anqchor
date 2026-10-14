# Welcome to qgenlib documentation

This is the documentation for the qgenlib project. 

## What is qgenlib?

qgenlib provides a set of classes and functions for setting up basic C++ programming environment for [Hyun Min Kang](https://github.com/hyunminkang). These libraries are designed with quantitative genetics in mind, but does not provide any sophisticated implementations of quantitative genetics algorithms. It is designed to be used with the [htslib](https://www.htslib.org/) library. Specifically, it provides the following functionalities:

* A macro-based command line parsing library implemented by Gonçalo Abecasis and used in many [tools](https://github.com/statgen).
* A collection of helper functions to output error, warning, or other messages.
* A collection of helper functions to handle HTS file formats supported by [htslib](https://www.htslib.org/).
* Classes and libraries to aid basic data wrangling using tab-limited files.
* Other utility functions to handle strings.
* A simple class to handle phred-scale quality scores.
* A simple class to handle genomic intervals.

## Documentation Overview

This documentation provides the following information:

* [Quickstart](quickstart.md): A quick guide to get started with qgenlib.
* [Install](install.md): How to install qgenlib.
* [Examples](examples.md): A simple example to demonstrate how to use qgenlib.