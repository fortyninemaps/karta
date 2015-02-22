import unittest
import numpy as np
from karta.vector.metadata import Metadata

class TestMetadata(unittest.TestCase):

    def setUp(self):
        self.onefield = Metadata(data=np.arange(200), singleton=False)
        self.multifield = Metadata(
                            data={"a":np.arange(200), "b":np.arange(200,400), 
                                  "c":np.arange(200)**2},
                            singleton=False)
        return

    def test_indexing(self):
        self.assertEqual(self.onefield[10], {"values": 10})
        self.assertEqual(self.multifield[10], {"a":10, "b": 210, "c": 100})
        return

    def test_slicing(self):
        self.assertTrue(np.all(self.onefield[5:15]["values"] == np.arange(5, 15)))
        res = self.multifield[5:15]
        self.assertTrue(np.all(res["a"] == np.arange(5, 15)))
        self.assertTrue(np.all(res["b"] == np.arange(205, 215)))
        self.assertTrue(np.all(res["c"] == np.arange(5, 15)**2))
        return

    def test_keys(self):
        self.assertTrue(np.all(self.onefield["values"] == np.arange(200)))
        self.assertTrue(np.all(self.multifield["a"] == np.arange(200)))
        self.assertTrue(np.all(self.multifield["b"] == np.arange(200, 400)))
        self.assertTrue(np.all(self.multifield["c"] == np.arange(200)**2))
        return

    def test_reference_vs_value_list(self):
        # metadata objects should carry values rather than references
        L = [1,2,3,4,5]
        md = Metadata(L, copy=True)
        L[3] = -99
        self.assertEqual(md["values"], [1,2,3,4,5])
        return

    def test_reference_vs_value_dict(self):
        # metadata objects should carry values rather than references
        D = {"A":[1,2,3,4,5],
             "B":[6,7,8,9,10]}
        md = Metadata(D, copy=True)
        D["A"][3] = -99
        self.assertEqual(md["A"], [1,2,3,4,5])
        return

    def test_nan_values(self):
        D = {"A":[1,np.nan,3,4,5],
             "B":[6,7,8,np.nan,10]}
        md = Metadata(D)
        self.assertTrue(np.isnan(md["A"][1]))
        self.assertEqual(md["A"], md["A"])
        return

    def test_none_values(self):
        D = {"A":[1,None,3,4,5],
             "B":[6,7,8,None,10]}
        md = Metadata(D)
        self.assertTrue(md["A"][1] is None)
        self.assertEqual(md["A"], md["A"])
        return

if __name__ == "__main__":
    unittest.main()

