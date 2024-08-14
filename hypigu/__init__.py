#
#   Copyright 2021--2024 Joshua Maglione
#
#   Distributed under MIT License
#

__version__ = 1.5

from .braid import BraidArrangementIgusa
from .constructors import CoxeterArrangement, LinialArrangement, ShiArrangement, CatalanArrangement, DirectSum, PolynomialToArrangement, ResonanceArrangement
from .gen_functions import FlagHilbertPoincareSeries, CoarseFlagHilbertPoincareSeries, IgusaZetaFunction
from .globals import verbose, ncpus
from .graded_poset import GradedPoset


__all__ = [
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
	'CoarseFlagHilbertPoincareSeries',
	'IgusaZetaFunction',
	'verbose',
	'ncpus',
].sort()
