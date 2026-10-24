"""Tests for debbie_images.py, offline, on a made-up catalogue.

    python3 -m unittest tools/test_debbie_images.py
"""
import contextlib
import csv
import gzip
import io
import math
import os
import random
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import debbie_images as di                                                   # noqa: E402

HEADER = ["name", "number", "ra", "dec", "part", "type", "shape_r", "mu_eff_g",
          "deepshadows_mean", "centroid_flag"]
JPEG = b"\xff\xd8\xff\xe0" + b"\0" * 32


def fake_catalogue(path, n=1000, seed=0):
    """n rows; every tenth has no score; the first two sit at dec 60."""
    rng = random.Random(seed)
    with gzip.open(path, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(HEADER)
        for k in range(n):
            num = 3 + 2 * k
            ra, dec = (150 + rng.uniform(-5, 5), rng.uniform(-5, 5))
            if k == 0:
                ra, dec = 100 + 0.9 / 60 / math.cos(math.radians(60)), 60.0   # 0.9' from (100, 60)
            if k == 1:
                ra, dec = 100 + 1.1 / 60 / math.cos(math.radians(60)), 60.0   # 1.1'
            score = "nan" if k % 10 == 9 else f"{(k % 100) / 100:.3f}"
            w.writerow([f"DEBBIE {num}", num, f"{ra:.7f}", f"{dec:.7f}",
                        ("original", "recovered", "extension")[k % 3], "SER",
                        f"{3 + (k % 50) / 5:.3f}", "25.0", score, 0])


def run(argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        di.main(argv)
    return out.getvalue()


class Selection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = os.path.join(cls.tmp.name, "cat.csv.gz")
        fake_catalogue(cls.path)
        cls.cat = di.Catalogue(cls.path)
        cls.all = list(range(len(cls.cat)))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_where_numeric_leaves_out_nan(self):
        got = di.where(self.cat, self.all, "deepshadows_mean>=0")
        self.assertEqual(len(got), 900)

    def test_where_text(self):
        self.assertEqual(len(di.where(self.cat, self.all, "part==extension")), 333)
        with self.assertRaises(SystemExit):
            di.where(self.cat, self.all, "part>extension")
        with self.assertRaises(SystemExit):
            di.where(self.cat, self.all, "no_such_column>1")

    def test_cone_uses_true_angles(self):
        got = di.cone(self.cat, self.all, 100.0, 60.0, 1.0)
        self.assertEqual(got, [0])

    def test_names(self):
        for n in ("DEBBIE 3", "DEBBIE_3", "3", "debbie 000003"):
            self.assertEqual(di.number_of(n), 3)
        with self.assertRaises(SystemExit):
            di.number_of("NGC 1")
        with contextlib.redirect_stdout(io.StringIO()) as out:
            got = di.by_number(self.cat, self.all, [3, 7, 4])
        self.assertEqual(got, [0, 2])
        self.assertIn("DEBBIE 4", out.getvalue())

    def test_stratified_even_and_weights_add_up(self):
        sel = di.where(self.cat, self.all, "deepshadows_mean>=0")
        with contextlib.redirect_stdout(io.StringIO()):
            got, cell, weight = di.stratified(self.cat, sel, ["deepshadows_mean:4"], 100, random.Random(1))
        self.assertEqual(len(got), 100)
        per = {}
        for i in got:
            per[cell[i]] = per.get(cell[i], 0) + 1
        self.assertEqual(sorted(per.values()), [25, 25, 25, 25])
        self.assertAlmostEqual(sum(weight.values()), len(sel))

    def test_stratified_small_cell_gives_its_share_to_others(self):
        sel = di.where(self.cat, self.all, "deepshadows_mean>=0")
        with contextlib.redirect_stdout(io.StringIO()):
            got, cell, weight = di.stratified(self.cat, sel, ["deepshadows_mean=0,0.03,1"], 100,
                                              random.Random(1))
        small = [i for i in got if cell[i] == "0<=deepshadows_mean<0.03"]
        self.assertEqual(len(got), 100)
        self.assertEqual(len(small), 30)          # scores 0, 0.01, 0.02: all 30 rows taken
        self.assertAlmostEqual(sum(weight.values()), len(sel))

    def test_sort_puts_nan_last(self):
        got = di.order(self.cat, self.all[:20], "deepshadows_mean", desc=True)
        col = self.cat.numeric("deepshadows_mean")
        vals = [col[i] for i in got]
        self.assertEqual(vals[:3], [0.18, 0.17, 0.16])
        self.assertTrue(math.isnan(vals[-1]) and math.isnan(vals[-2]))
        args = di.parser().parse_args(["--sort", "shape_r", "--desc"])
        self.assertTrue(args.desc)

    def test_field(self):
        a = di.parser().parse_args([])
        self.assertEqual(di.field(3.0, a), (153, 0.262))
        self.assertEqual(di.field(18.0, a), (550, 0.262))
        a.fov_reff = 1000
        self.assertEqual(di.field(18.0, a)[0], 3000)

    def test_valid(self):
        self.assertTrue(di.valid(JPEG, "jpg"))
        self.assertTrue(di.valid(b"SIMPLE  =                    T", "fits"))
        self.assertFalse(di.valid(b"no such layer", "jpg"))
        self.assertFalse(di.valid(b"<html>busy</html>", "fits"))


class Download(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "cat.csv.gz")
        fake_catalogue(self.path, n=60)
        self.out = os.path.join(self.tmp.name, "out")

    def tearDown(self):
        self.tmp.cleanup()

    def answer(self, body):
        resp = mock.MagicMock()
        resp.__enter__.return_value.read.return_value = body
        return mock.patch.object(di.urllib.request, "urlopen", return_value=resp)

    def test_fetch_refuses_a_text_answer(self):
        dest = os.path.join(self.tmp.name, "x.jpg")
        with self.answer(b"no such layer"), mock.patch.object(di.time, "sleep"):
            why = di.fetch("http://x", dest, "jpg")
        self.assertIn("not a jpg", why)
        self.assertFalse(os.path.exists(dest))

    def test_run_writes_images_gallery_and_resumes(self):
        argv = ["--catalogue", self.path, "--out", self.out, "--sample", "12",
                "--stratify", "shape_r:3", "--context"]
        with self.answer(JPEG), mock.patch.object(di.time, "sleep"):
            first = run(argv)
            again = run(argv)
        self.assertIn(", 48 cutouts to fetch", first)          # 12 entries x 4 panels
        self.assertIn(", 0 cutouts to fetch", again)
        rows = list(csv.DictReader(open(os.path.join(self.out, "selection.csv"))))
        self.assertEqual(len(rows), 12)
        page = open(os.path.join(self.out, "index.html")).read()
        for r in rows:
            d = os.path.join(self.out, r["name"].replace(" ", "_"))
            self.assertEqual(sorted(os.listdir(d)), ["context.jpg", "image.jpg", "model.jpg", "resid.jpg"])
            self.assertIn(r["name"], page)

    def test_cap(self):
        with self.assertRaises(SystemExit):
            run(["--catalogue", self.path, "--out", self.out, "--max-objects", "10"])


if __name__ == "__main__":
    unittest.main()
