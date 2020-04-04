# RODEO2antiSMASH
## Overview
This is a frozen version of a script used in our [paper](https://www.frontiersin.org/articles/10.3389/fgene.2020.00226/full). If you want to use it for your project, better consider [the python package](https://github.com/bikdm12/rodeo_utils). It contains all functionality of this script as well as some other tools.

This notebook modifies GenBank files to make them look like [antiSMASH](https://antismash.secondarymetabolites.org/) results. The additional information is taken from the output of [RODEO](http://ripp.rodeo/). The resulting GenBank files can be used to construct a sequence similarity network with [BiG-SCAPE](https://bigscape-corason.secondarymetabolites.org/)
## Requirenments
* Python 2.7
* Biopython
