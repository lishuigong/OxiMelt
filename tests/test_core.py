import unittest
from oximelt.core import predict_formula, OxiMeltError, parse_formula

class TestOxiMelt(unittest.TestCase):
    def assertClose(self, got, expected, tol=0.05):
        self.assertLessEqual(abs(got - expected), tol)

    def test_luag_global(self):
        r = predict_formula("Lu3Al5O12")
        self.assertClose(r.predicted_Tm_K, 2309.75586)
        self.assertFalse(r.hf_local_correction_applied)
        self.assertEqual(r.prediction_mode, "GLOBAL 3P")

    def test_yb2tio5_global(self):
        self.assertClose(predict_formula("Yb2TiO5").predicted_Tm_K, 2336.52918)

    def test_yb2ti2o7_global(self):
        self.assertClose(predict_formula("Yb2Ti2O7").predicted_Tm_K, 2220.72033)

    def test_parentheses(self):
        c = parse_formula("(Lu0.5Y0.5)2O3")
        self.assertAlmostEqual(c["Lu"], 1.0)
        self.assertAlmostEqual(c["Y"], 1.0)
        self.assertAlmostEqual(c["O"], 3.0)

    def test_peroxide_rejected(self):
        with self.assertRaises(OxiMeltError):
            predict_formula("Na2O2")


    def test_hf_lu_binary_oof_calibration(self):
        # Exact frozen Lu support exists at x_Hf = 0.60.
        r = predict_formula("Hf0.6Lu0.4O1.8")
        self.assertClose(r.global_Tm_K, 2739.317405, tol=0.06)
        self.assertClose(r.local_correction_K, 456.585923, tol=0.06)
        self.assertClose(r.predicted_Tm_K, 3195.903328, tol=0.08)
        self.assertTrue(r.hf_local_domain_eligible)
        self.assertTrue(r.hf_local_correction_applied)
        self.assertAlmostEqual(r.hf_local_support_coverage, 1.0)
        self.assertIn("OOF CALIBRATED", r.applicability)
        self.assertIn("Hf-BINARY OOF CALIBRATION", r.prediction_mode)

    def test_hf_lu_sc_75_25(self):
        r = predict_formula("Hf0.6Lu0.3Sc0.1O1.8")
        self.assertClose(r.global_Tm_K, 2713.726975)
        self.assertClose(r.local_correction_K, 433.691361)
        self.assertClose(r.predicted_Tm_K, 3147.418336)
        self.assertTrue(r.hf_local_correction_applied)
        self.assertAlmostEqual(r.hf_local_support_coverage, 1.0)

    def test_hf_lu_sc_50_50(self):
        r = predict_formula("Hf0.6Lu0.2Sc0.2O1.8")
        self.assertClose(r.global_Tm_K, 2703.806313)
        self.assertClose(r.local_correction_K, 410.796799)
        self.assertClose(r.predicted_Tm_K, 3114.603112)
        self.assertTrue(r.hf_local_correction_applied)

    def test_hf_lu_sc_25_75(self):
        r = predict_formula("Hf0.6Lu0.1Sc0.3O1.8")
        self.assertClose(r.global_Tm_K, 2709.554988)
        self.assertClose(r.local_correction_K, 387.902237)
        self.assertClose(r.predicted_Tm_K, 3097.457225)
        self.assertTrue(r.hf_local_correction_applied)

    def test_hf_lu_yb_75_25(self):
        r = predict_formula("Lu0.3Yb0.1Hf0.6O1.8")
        self.assertClose(r.local_correction_K, 452.260066)
        self.assertClose(r.predicted_Tm_K, 3170.751432)
        self.assertTrue(r.hf_local_correction_applied)

    def test_hf_lu_yb_50_50(self):
        r = predict_formula("Lu0.2Yb0.2Hf0.6O1.8")
        self.assertClose(r.local_correction_K, 447.934209)
        self.assertClose(r.predicted_Tm_K, 3161.224476)
        self.assertTrue(r.hf_local_correction_applied)

    def test_hf_lu_yb_25_75(self):
        r = predict_formula("Lu0.1Yb0.3Hf0.6O1.8")
        self.assertClose(r.local_correction_K, 443.608352)
        self.assertClose(r.predicted_Tm_K, 3167.322028)
        self.assertTrue(r.hf_local_correction_applied)

    def test_hf_family_insufficient_support_falls_back_to_global(self):
        # Y has no frozen binary residual profile in the published local library.
        r = predict_formula("Hf0.6Lu0.1Y0.3O1.8")
        self.assertTrue(r.hf_local_domain_eligible)
        self.assertFalse(r.hf_local_correction_applied)
        self.assertLess(r.hf_local_support_coverage, 0.80)
        self.assertClose(r.predicted_Tm_K, r.global_Tm_K, tol=1e-9)

if __name__ == "__main__":
    unittest.main()
