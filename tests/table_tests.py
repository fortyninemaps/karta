import unittest
import json
import numpy as np
from karta.vector.table import Table, merge

class TestTable(unittest.TestCase):

    def setUp(self):
        self.onefield = Table(data=np.arange(200))
        self.multifield = Table(data={"a":np.arange(200),
                                         "b":np.arange(200,400),
                                         "c":np.arange(200)**2})

    def test_asjson_onefield(self):
        d = json.loads(self.onefield.asjson())
        self.assertEqual(len(d), 1)
        self.assertTrue("value" in d)

    def test_asjson_multifield(self):
        d = json.loads(self.multifield.asjson())
        self.assertEqual(len(d), 3)
        self.assertTrue("a" in d)
        self.assertTrue("b" in d)
        self.assertTrue("c" in d)

    def test_indexing(self):
        self.assertEqual(self.onefield[10], (10,))
        self.assertEqual(set(self.multifield[10]), set([10, 210, 100]))

    def test_get(self):
        self.assertEqual(self.onefield.get(10), {"value": 10})
        self.assertEqual(self.multifield.get(10), {"a":10, "b": 210, "c": 100})

    def test_slicing(self):
        self.assertTrue(np.all(self.onefield.get(slice(5,15))["value"] == np.arange(5, 15)))
        res = self.multifield.get(slice(5,15))
        self.assertTrue(np.all(res["a"] == np.arange(5, 15)))
        self.assertTrue(np.all(res["b"] == np.arange(205, 215)))
        self.assertTrue(np.all(res["c"] == np.arange(5, 15)**2))

    def test_keys(self):
        self.assertTrue(np.all(self.onefield.getfield("value") == np.arange(200)))
        self.assertTrue(np.all(self.multifield.getfield("a") == np.arange(200)))
        self.assertTrue(np.all(self.multifield.getfield("b") == np.arange(200, 400)))
        self.assertTrue(np.all(self.multifield.getfield("c") == np.arange(200)**2))

    def test_reference_vs_value_list(self):
        # metadata objects should carry values rather than references
        L = [1,2,3,4,5]
        md = Table(L)
        L[3] = -99
        self.assertEqual(md.getfield("value"), [1,2,3,4,5])

    def test_reference_vs_value_dict(self):
        # metadata objects should carry values rather than references
        D = {"A":[1,2,3,4,5],
             "B":[6,7,8,9,10]}
        md = Table(D)
        D["A"][3] = -99
        self.assertEqual(md.getfield("A"), [1,2,3,4,5])

    def test_nan_values(self):
        D = {"A":[1,np.nan,3,4,5],
             "B":[6,7,8,np.nan,10]}
        md = Table(D)
        self.assertTrue(np.isnan(md.getfield("A")[1]))
        self.assertEqual(md.getfield("A"), md.getfield("A"))

    def test_none_values(self):
        D = {"A":[1,None,3,4,5],
             "B":[6,7,8,None,10]}
        md = Table(D)
        self.assertTrue(md.getfield("A")[1] is None)
        self.assertEqual(md.getfield("A"), md.getfield("A"))

    def test_init_specified_fields(self):
        md = Table([(1,2,3),(4,5,6),(7,8,9)], fields=("a","b","c"))
        self.assertEqual(md.data, [(1,2,3),(4,5,6),(7,8,9)])
        self.assertEqual(md.fields, set(("a","b","c")))

    def test_init_specified_fields2(self):
        md = Table([("by air",),("by land,"),("by sea",)], fields=("mode",))
        self.assertEqual(md.data, [("by air",),("by land,"),("by sea",)])
        self.assertEqual(md.fields, set(("mode",)))

    def test_extend(self):
        md1 = Table({"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                        "number": [3784, 78, 83]})
        md2 = Table({"street": ["Winding Ln.", "Oyster Pond Rd."],
                        "color": ["Grey", "White"]})
        md1.extend(md2)
        self.assertEqual(len(md1), 5)
        self.assertEqual(set((None, "Winding Ln.")), set(md1[3]))

    def test_addfield(self):
        # data, field syntax used in order to specify field order
        md = Table([("F R. Lillie", 3784), ("School St.", 78), ("Quissett Ave.", 83)],
                      fields=("street", "number"))
        md.setfield("driveway", [True, True, False])
        ans = Table([("F R. Lillie", 3784, True), ("School St.", 78, True),
                        ("Quissett Ave.", 83, False)],
                       fields=("street", "number", "driveway"))
        self.assertEqual(md, ans)

    def test_updatefield(self):
        md = Table({"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                       "number": [3784, 78, 83]})
        md.setfield("number", [3782, 78, 83])
        ans = Table({"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                        "number": [3782, 78, 83]})
        self.assertEqual(md, ans)

    def test_update_from_table(self):
        tbl1 = Table({"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                       "number": [3784, 78, 83]})
        tbl2 = Table({"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                       "number": [3785, 79, 82],
                       "color": ["red", "green", "orange"]})

        new = tbl1.updated(tbl2)

        self.assertTrue("color" in new.fields)
        self.assertEqual(new.getfield("number"), [3785, 79, 82])

        self.assertFalse("color" in tbl1.fields)
        self.assertEqual(tbl1.getfield("number"), [3784, 78, 83])

    def test_update_from_dict(self):
        tbl1 = Table({"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                       "number": [3784, 78, 83]})
        tbl2 = {"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                "number": [3785, 79, 82],
                "color": ["red", "green", "orange"]}

        new = tbl1.updated(tbl2)

        self.assertTrue("color" in new.fields)
        self.assertEqual(new.getfield("number"), [3785, 79, 82])

        self.assertFalse("color" in tbl1.fields)
        self.assertEqual(tbl1.getfield("number"), [3784, 78, 83])

    def test_merge(self):
        tbl1 = Table({"street": ["F. R. Lillie", "School St.", "Quissett Ave."],
                      "number": [3784, 78, 83]})
        tbl2 = Table({"street": ["Winding Ln.", "Oyster Pond Rd."],
                      "number": [17, 33]})
        d1 = {"street": "Fern Lane", "number": 33}

        merged = merge([tbl1, d1, tbl2])
        self.assertEqual(len(merged), 6)
        self.assertEqual(merged.fields, set(("street", "number")))
        self.assertEqual(set(merged[3]), set(("Fern Lane", 33)))


if __name__ == "__main__":
    unittest.main()

