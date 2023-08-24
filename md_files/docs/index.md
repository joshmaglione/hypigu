# HypIgu 

Author: [Joshua Maglione](https://www.math.uni-bielefeld.de/~jmaglione/).

Documentation for the [HypIgu](https://github.com/joshmaglione/hypigu) package for [SageMath](https://www.sagemath.org/).

## Purpose

The goal of **HypIgu** (**HYP**erplane **IGU**sa) is to provide SageMath with the functionality to compute various zeta functions associated with hyperplane arrangements. Included are common constructions for hyperplane arrangements and specializations of the flag Hilbert&ndash;Poincar&#233; series defined in [Maglione&ndash;Voll](https://arxiv.org/abs/2103.03640). Mathematical details are given in [Maglione&ndash;Voll](https://arxiv.org/abs/2103.03640). We outline the functions included in HypIgu and provide example cases. 

## Setup

The simplest way to **install** HypIgu is to run the following 

```sh
$ sage --pip install hypigu
```

Alternatively, one can download the [latest release](https://github.com/joshmaglione/hypigu/releases/latest) and unzip it into a directory that SageMath can find for importing.

To **update** an older version of HypIgu to the latest version, run the following 

```sh 
$ sage --pip install hypigu --upgrade 
```

HypIgu has no external dependencies and is compatible with SageMath 9.6 and later. It may work just fine with earlier versions of SageMath, but these have not been tested.

## Importing

Import HypIgu during your SageMath run with the following

```python
import hypigu as hi
```

Throughout this documentation, we use `hi` for the reference name of `hypigu`.

## Funding 

This work was supported in part by DFG-grant [373111162](https://gepris.dfg.de/gepris/projekt/373111162?language=en).

## References 

1. [Joshua Maglione](https://www.math.uni-bielefeld.de/~jmaglione/) and [Christopher Voll](https://www.math.uni-bielefeld.de/~voll/). Flag Hilbert&ndash;Poincar&#233; series of hyperplane arrangements and Igusa zeta functions, to appear in Israel J. Math. 2021. [arXiv:2103:03640](https://arxiv.org/abs/2103.03640).
   