import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch, MagicMock

from oximelt import core


# Synthetic formulas below test the calculation, not experimental observations.
class TestVacancyScaling(unittest.TestCase):
    def test_synthetic_nominal_stoichiometry(self):
        r = core.predict_formula('Hf0.9Sc0.1O1.95')
        s = r.local_support_details[0]
        self.assertAlmostEqual(r.nominal_oxygen_vacancy_site_fraction, 0.025)
        self.assertAlmostEqual(s.support_vacancy_fraction, .1)
        self.assertAlmostEqual(s.vacancy_scaling_factor, .25)
        self.assertAlmostEqual(s.dopant_share, 1)
        self.assertAlmostEqual(r.local_correction_K, .25 * 365.0076754188385)
        self.assertAlmostEqual(r.predicted_Tm_C, r.global_Tm_C + .25 * 365.0076754188385)
        self.assertTrue(r.hf_local_concentration_extrapolation)
        self.assertEqual(s.mode, 'vacancy_scaled_extrapolation')

    def test_formula_scale_invariance(self):
        a = core.predict_formula('Hf0.9Sc0.1O1.95')
        b = core.predict_formula('Hf18Sc2O39')
        self.assertAlmostEqual(a.global_Tm_K, b.global_Tm_K)
        self.assertAlmostEqual(a.local_correction_K, b.local_correction_K)

    def test_unicode_input_same_result(self):
        a = core.predict_formula('Hf₀.₉Sc₀.₁O₁.₉₅')
        self.assertAlmostEqual(a.predicted_Tm_C, core.predict_formula('Hf0.9Sc0.1O1.95').predicted_Tm_C)

    def test_sc_half_dopant_fraction_halves_residual(self):
        a = core.predict_formula('Hf0.9Sc0.1O1.95')
        b = core.predict_formula('Hf0.95Sc0.05O1.975')
        self.assertAlmostEqual(a.local_correction_K, 2*b.local_correction_K)

    def test_small_positive_vacancy(self):
        r = core.predict_formula('Hf0.996Sc0.004O1.998')
        self.assertAlmostEqual(r.nominal_oxygen_vacancy_site_fraction, .001)
        self.assertTrue(r.hf_local_correction_applied)
        self.assertAlmostEqual(r.local_correction_K, .01*365.0076754188385)

    def test_lu_uses_nearest_vacancy_support_not_interpolation(self):
        r = core.predict_formula('Hf0.76Lu0.24O1.88')
        s = r.local_support_details[0]
        self.assertEqual(s.support_formula, '0.15Lu2O3-0.85HfO2')
        self.assertAlmostEqual(r.local_correction_K, .06/(.15/2.3)*362.8404862994107)
        self.assertEqual(s.mode, 'vacancy_scaled_nearest')
        self.assertFalse(s.concentration_extrapolation)

    def test_equal_distance_selects_lower_vacancy_fraction(self):
        rows = [dict(dopant='Sc', x_Hf=1-4*f, residual_K=res,
                     n_records_at_x=1, support_formula=name,
                     dopant_formal_valence=3, nominal_vacancy_site_fraction=f)
                for f, res, name in [(.1,900,'upper'),(.08,100,'lower')]]
        with patch.dict(core.HF_LOCAL_SUPPORT, {'Sc': rows}):
            r = core.predict_formula('Hf0.64Sc0.36O1.82')
        self.assertEqual(r.local_support_details[0].support_formula, 'lower')
        self.assertAlmostEqual(r.local_correction_K, .09/.08*100)

    def test_duplicate_support_fraction_uses_median(self):
        base = dict(core.HF_LOCAL_SUPPORT['Sc'][0])
        rows = [dict(base, residual_K=v, support_formula=str(v)) for v in [900,100,200]]
        with patch.dict(core.HF_LOCAL_SUPPORT, {'Sc': rows}):
            r = core.predict_formula('Hf0.9Sc0.1O1.95')
        self.assertAlmostEqual(r.local_support_details[0].unscaled_residual_K, 200)
        self.assertAlmostEqual(r.local_correction_K, .25*200)
        self.assertEqual(r.local_support_details[0].support_n, 3)

    def test_signed_residual_is_not_clipped(self):
        row = dict(core.HF_LOCAL_SUPPORT['Sc'][0], residual_K=-100)
        with patch.dict(core.HF_LOCAL_SUPPORT, {'Sc': [row]}):
            r = core.predict_formula('Hf0.9Sc0.1O1.95')
        self.assertAlmostEqual(r.local_correction_K, -25.0)

    def test_valence_mismatch_is_not_transferred(self):
        row = dict(core.HF_LOCAL_SUPPORT['Sc'][0], dopant_formal_valence=2)
        with patch.dict(core.HF_LOCAL_SUPPORT, {'Sc': [row]}):
            r = core.predict_formula('Hf0.9Sc0.1O1.95')
        self.assertFalse(r.hf_local_correction_applied)

    def test_coverage_boundary_without_renormalization(self):
        r = core.predict_formula('Hf0.6Sc0.32Y0.08O1.8')
        self.assertAlmostEqual(r.hf_local_support_coverage, .8)
        self.assertTrue(r.hf_local_correction_applied)
        self.assertAlmostEqual(r.local_correction_K, .8*365.0076754188385)
        r = core.predict_formula('Hf0.6Sc0.3196Y0.0804O1.8')
        self.assertFalse(r.hf_local_correction_applied)
        self.assertEqual(r.local_correction_K, 0)

    def test_outside_family_remains_global(self):
        for f in ['HfO2', 'Hf0.6Zr0.4O2', 'Hf0.54Sc0.46O1.77',
                  'K8Hf(MoO4)6', 'K2Hf(MoO4)3', 'Hf0.6Lu0.1Sc0.1Yb0.2O1.8']:
            with self.subTest(formula=f):
                r = core.predict_formula(f)
                self.assertFalse(r.hf_local_correction_applied)
                self.assertEqual(r.predicted_Tm_K, r.global_Tm_K)

    def test_mixed_batch_exports_all_audit_fields(self):
        for formulas in [['bad','Hf0.9Sc0.1O1.95'],['Hf0.9Sc0.1O1.95','bad']]:
            with TemporaryDirectory() as td:
                p = Path(td)
                (p/'input.csv').write_text('formula\n'+'\n'.join(formulas), encoding='utf8')
                self.assertEqual(core.predict_csv(p/'input.csv',p/'out.csv'),2)
                with (p/'out.csv').open(encoding='utf-8-sig',newline='') as f:
                    rows = list(csv.DictReader(f))
                good = next(r for r in rows if not r['error'])
                support = json.loads(good['hf_local_support_details_json'])[0]
                self.assertAlmostEqual(support['vacancy_scaling_factor'], .25)
                self.assertEqual(good['hf_local_concentration_extrapolation'],'1')

    def test_gui_prediction_handler_reports_scaling(self):
        try:
            from oximelt.gui import OxiMeltApp
        except ImportError:
            self.skipTest('Tkinter is not available on this Python runtime')
        # Exercise the GUI callback without opening a desktop window.
        app = MagicMock()
        app.formula_var.get.return_value = 'Hf0.9Sc0.1O1.95'
        OxiMeltApp._predict_single(app)
        summary = app._set_summary.call_args.args[0]
        self.assertIn('ratio=0.2500',summary)
        self.assertIn('Concentration extrapolation: Yes',summary)
        self.assertIn('Family estimate',app.result_c.config.call_args.kwargs['text'])


if __name__ == '__main__':
    unittest.main()
