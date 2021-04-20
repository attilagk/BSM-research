# BSM research

## What's this?

* my computational research project on brain somatic mosaicism (BSM)
* the history of related statistical DNA analysis a.k.a my **lab notebook**
* collection of scripts and small data sets for project-related computation
* slides for scientific presentations

## Introduction

The present research project has been conducted at the Mount Sinai School of Medicine in Andy Chess' laboratory, which is part of an NIH funded research consortium called Brain Somatic Mosaicism Network (BSMN).  Short, general information on the BSMN is available on its website [http://brainsomaticmosaicism.org/](http://brainsomaticmosaicism.org/).  A more extensive document of general information is the review article authored by the BSMN ([Science. 2017 Apr 28;356(6336). pii: eaal1641.](http://science.sciencemag.org/content/356/6336/eaal1641)).

Within the BSMN this research project is specific in that it investigates the role of somatic mosaicism in schizophrenia somatic variants are called from deep DNA sequencing data from sample pairs.  The paired sample approach, originally developed in cancer research, allows a more accurate calling of somatic variants in the face of germline variants.

## Structure and content

```
ndar/           # for data sharing at https://ndar.nih.gov/
notebook/       # the main matter: my lab notebook
presentations/  # LaTeX/Beamer slides
src/            # scripts for project-wide use
tables/         # small data sets compiled by hand
```

As noted above the main matter of this repository is the lab notebook, a chronological sequence of articles in the `notebook/` directory.  Technically these articles are R Markdown documents or iPython notebooks, which means that they contain code chunks to be executed by interpreters such as `R`, `python` or `bash`.  Longer, more complex code chunks have been moved to scripts located either in `notebook/` (these are associated to some article) or in `src/` (these are shared by multiple articles).  Finally, the most general-purpose code is part of [research-tools](https://github.com/attilagk/research-tools), a separate repository.

## The data

Of course, the data---Terabytes of DNA sequences---are not contained in this repository but they do exist on our computational server and at the [NIMH Data Archive](https://ndar.nih.gov/) (NDA), which will release them for public use.

## Reproduce

* clone this git project `https://github.com/attilagk/BSM-research`
* clone `https://github.com/attilagk/src`
* modules needed
    * argparse
    * copy
    * ensembl_rest
    * functools
    * io
    * itertools
    * numpy
    * operator
    * os.path
    * pandas
    * pickle
    * re
    * scipy.stats
    * statsmodels
    * subprocess
    * synapseclient
* set `$PYTHONPATH` in your `.profile` or `.bashrc` as follows
```
export PYTHONPATH="/path/to/bsm/src:$PYTHONPATH"
```
