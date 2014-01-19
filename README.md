comcluster
==========

Selects a representative set for a publishing community (e.g. editorial board for a journal).

### Setup

Beforing using ```cluster.py``` install [graph-tool](http://graph-tool.skewed.de/).

### Usage

Use the documentation included in the python script by running

``` python cluster.py -h ```

### Example - Principles of Programming Language (POPL) Conference

The supplied ```popl-sample/``` contains:

- ```data/``` containing sample POPL data generated using [pubminer](https://github.com/rjullman/pubminer) that can be used as input to ```miner.py```;

- ```output/``` containing sample text and graphical output generated from the POPL community; and

- ```actual/``` containing the actual program committees for various past years of POPL.  This serves as an interesting comparison point for generated results. 

To get started try running

```
python cluster.py popl-sample/data/popl2013.dat -a 4 -n 30 -w 3
```
