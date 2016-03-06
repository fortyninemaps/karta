import unittest
from karta.raster import _dem

class DEMDriverTests(unittest.TestCase):

    def test_nreps_re(self):
        nr = _dem.nreps_re("2(I4,I2,F7.4)")
        self.assertEqual(nr, 2)
        nr = _dem.nreps_re("2(3I6)")
        self.assertEqual(nr, 6)
        nr = _dem.nreps_re("4(6F3.7)")
        self.assertEqual(nr, 24)
        return

    def test_parse(self):
        p = _dem.parse("2(I4,I2,F7.4)", ' -81 0 0.0000  82 0 0.0000')
        self.assertEqual(p, [-81, 0, 0.0, 82, 0, 0.0])
        return

    def test_dtype(self):
        t = _dem.dtype("2(I4,D24.15,F7.4)")
        self.assertEqual(t, [int, _dem.coerce_float, _dem.coerce_float])
        return

    def test_reclen(self):
        n = _dem.reclen("2(I4,I2,F7.4)")
        self.assertEqual(n, [4, 2, 7])
        return

if __name__ == "__main__":
    unittest.main()
