import unittest
from hypigu.constructors import _Coxeter_check, CoxeterArrangement  

class TestHyperplaneArrangements(unittest.TestCase):

	def test_Coxeter_check_valid(self):
		self.assertTrue(all(_Coxeter_check('A', k) for k in range(1, 101)))
		self.assertTrue(all(_Coxeter_check('B', k) for k in range(1, 101)))
		self.assertTrue(all(_Coxeter_check('C', k) for k in range(1, 101)))
		self.assertTrue(all(_Coxeter_check('D', k) for k in range(1, 101)))
		self.assertTrue(_Coxeter_check('E', 8))
		self.assertTrue(_Coxeter_check('G', 2))
		self.assertTrue(_Coxeter_check('H', 4))
		self.assertTrue(_Coxeter_check('F', 4))
		self.assertTrue(_Coxeter_check('I', 42))
		self.assertFalse(_Coxeter_check('A', -1))
		self.assertFalse(_Coxeter_check('A', 0))
		self.assertFalse(_Coxeter_check('Z', 3))
		self.assertFalse(_Coxeter_check('E', 9))
		self.assertFalse(_Coxeter_check('F', 5))
		self.assertFalse(_Coxeter_check('G', 3))
		self.assertFalse(_Coxeter_check('H', 1))

	# These are "loose" checks.
	def test_CoxeterArrangement(self):
		for m in range(3, 10):
			I2_m = CoxeterArrangement('I', m)
			L = I2_m.matroid().lattice_of_flats()
			self.assertEqual(L.rank(), 2)
			self.assertEqual(len(L), m + 2)
		for m in range(3, 5):
			for n in range(10, 12):
				A = CoxeterArrangement(f"I{m} I{n}")
				L = A.matroid().lattice_of_flats()
				self.assertEqual(L.rank(), 4)
				self.assertEqual(len(L.atoms()), m + n)
		for m in range(1, 5):
			for n in range(2, 5):
				A = CoxeterArrangement(f"A{m} B{n}")
				L = A.matroid().lattice_of_flats()
				self.assertEqual(L.rank(), m + n)
				self.assertEqual(len(L.atoms()), m*(m + 1)//2 + n*n)
		for m in range(1, 5):
			A = CoxeterArrangement(f"D{m}")
			L = A.matroid().lattice_of_flats()
			self.assertEqual(L.rank(), m)
			self.assertEqual(len(L.atoms()), m*(m - 1))
		

if __name__ == '__main__':
	unittest.main()
