import unittest
import fwdpy11.genetic_values
import numpy as np

# TODO reimplement to test new API

# class testSlocusAdditive(unittest.TestCase):
#     """
#     API and behavior tests
#     """
# 
#     def testConstructionWithNumpyFloat(self):
#         x = fwdpy11.genetic_values.SlocusAdditive(np.float(1.0))
#         self.assertEqual(x.scaling, 1.0)
#         self.assertEqual(x.is_fitness, True)
#         self.assertTrue(isinstance(x.gvalue_to_fitness,
#                                    fwdpy11.genetic_values.GeneticValueIsFitness))
#         self.assertTrue(isinstance(x.noise,
#                                    fwdpy11.genetic_value_noise.SlocusPopNoNoise))
# 
#     def testConstructWithNaN(self):
#         with self.assertRaises(ValueError):
#             x = fwdpy11.genetic_values.SlocusAdditive(np.nan)
# 
#     def testConstructTrait(self):
#         x = fwdpy11.genetic_values.SlocusAdditive(
#             1.0, fwdpy11.genetic_values.GSS(1.0, 0.0))
#         self.assertEqual(x.scaling, 1.0)
#         self.assertEqual(x.is_fitness, False)
#         self.assertTrue(isinstance(x.gvalue_to_fitness,
#                                    fwdpy11.genetic_values.GSS))
# 
# 
# class testSlocusMult(unittest.TestCase):
#     """
#     API and behavior tests
#     """
# 
#     def testConstructionWithNumpyFloat(self):
#         x = fwdpy11.genetic_values.SlocusMult(np.float(1.0))
#         self.assertEqual(x.scaling, 1.0)
#         self.assertEqual(x.is_fitness, True)
#         self.assertTrue(isinstance(
#             x.noise, fwdpy11.genetic_value_noise.SlocusPopNoNoise))
# 
#     def testConstructWithNaN(self):
#         with self.assertRaises(ValueError):
#             x = fwdpy11.genetic_values.SlocusMult(np.nan)
# 
#     def testConstructTrait(self):
#         x = fwdpy11.genetic_values.SlocusMult(
#             1.0, fwdpy11.genetic_values.GSS(1.0, 0.0))
#         self.assertEqual(x.scaling, 1.0)
#         self.assertEqual(x.is_fitness, False)
#         self.assertTrue(isinstance(x.gvalue_to_fitness,
#                                    fwdpy11.genetic_values.GSS))
# 
# 
# class testGBR(unittest.TestCase):
#     def testConstruct(self):
#         x = fwdpy11.genetic_values.GBR(fwdpy11.genetic_values.GSS(1.0, 0))
#         self.assertEqual(x.gvalue_to_fitness.VS, 1.0)
#         self.assertEqual(x.gvalue_to_fitness.opt, 0.0)
#         self.assertTrue(isinstance(
#             x.noise, fwdpy11.genetic_value_noise.SlocusPopNoNoise))
# 
#     def testBadConstruct(self):
#         with self.assertRaises(TypeError):
#             x = fwdpy11.genetic_values.GBR(
#                 fwdpy11.genetic_values.GeneticValueIsFitness())
# 
# 
# if __name__ == "__main__":
#     pass
