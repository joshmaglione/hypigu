# HypIgu 

Author: Joshua Maglione.

Documentation for the [HypIgu](https://github.com/joshmaglione/hypigu) package for [SageMath](https://www.sagemath.org/).

## Purpose

The goal of HypIgu (HYPerplane IGUsa) is to provide SageMath with the functionality to compute various zeta functions associated with hyperplane arrangements. Included are common constructions for hyperplane arrangements and specializations of the flag Hilbert&ndash;Poincar&#233; series defined in Maglione&ndash;Voll. Mathematical details are given in Maglione&ndash;Voll. We outline the functions included in HypIgu and provide example cases. 

## Setup

The simplest way to install HypIgu is to run the following 

```bash
$ sage --pip install hypigu
```

Alternatively, one can download the latest release and unzip it into a directory that SageMath can find for importing.

HypIgu has no external dependencies and is compatible with SageMath 9.2.

## Importing

Import HypIgu during your SageMath run with the following

```python
import hypigu as hi
```

Throughout we use `hi` for the reference name of `hypigu`.

## Funding 

This work was supported in part by DFG-grant [373111162](https://gepris.dfg.de/gepris/projekt/373111162?language=en).
