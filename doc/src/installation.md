# Installation

```
Pkg.add("PGFPlots")
```

In addition, you will need to install the following dependencies if you do not already have them on your system.
* Pdf2svg. This is required by TikzPictures. On Ubuntu, you can get this by running `sudo apt-get install pdf2svg` and on RHEL/Fedora by running `sudo dnf install pdf2svg`. On Windows, you can download the binaries from http://www.cityinthesky.co.uk/opensource/pdf2svg/. Be sure to add pdf2svg to your path (and restart).
* Pgfplots (version 1.10 or later). Install using your latex package manager (e.g., texlive or miktex).

Once these things are installed, you should be able to run the following:

```@example runningExample
using PGFPlots
```