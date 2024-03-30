#
#   Copyright 2021--2023 Joshua Maglione
#
#   Distributed under MIT License
#

__version__ = 1.5

from .braid import BraidArrangementIgusa
from .constructors import CoxeterArrangement, LinialArrangement, ShiArrangement, CatalanArrangement, DirectSum, PolynomialToArrangement, ResonanceArrangement
from .lattice_flats import LatticeOfFlats
from .gen_functions import FlagHilbertPoincareSeries, IgusaZetaFunction, CoarseFlagHPSeries, AnalyticZetaFunction, AtomZetaFunction, TopologicalZetaFunction
from .database import internal_database

__all__ = [
	'BraidArrangementIgusa',
	'CoxeterArrangement',
	'LinialArrangement',
	'ShiArrangement',
	'CatalanArrangement',
	'DirectSum',
	'PolynomialToArrangement',
	'ResonanceArrangement',
	'LatticeOfFlats',
	'FlagHilbertPoincareSeries',
	'IgusaZetaFunction',
	'CoarseFlagHPSeries',
	'AnalyticZetaFunction',
	'AtomZetaFunction',
	'TopologicalZetaFunction',
	'internal_database'
]
