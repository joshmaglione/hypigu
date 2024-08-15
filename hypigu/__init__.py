#
#   Copyright 2021--2024 Joshua Maglione
#
#   Distributed under MIT License
#

__version__ = 1.5

from .braid import BraidArrangementIgusa
from .constructors import CoxeterArrangement, LinialArrangement, ShiArrangement, CatalanArrangement, DirectSum, PolynomialToArrangement, ResonanceArrangement
from .gen_functions import FlagHilbertPoincareSeries, CoarseFHPSeries, IgusaZetaFunction, TopologicalZetaFunction
from .globals import verbose
from .graded_poset import GradedPoset


__all__ = sorted([
	'BraidArrangementIgusa',
	'CoxeterArrangement',
	'LinialArrangement',
	'ShiArrangement',
	'CatalanArrangement',
	'DirectSum',
	'PolynomialToArrangement',
	'ResonanceArrangement',
	'GradedPoset',
	'FlagHilbertPoincareSeries',
	'CoarseFHPSeries',
	'IgusaZetaFunction',
	'verbose',
	'TopologicalZetaFunction',
])