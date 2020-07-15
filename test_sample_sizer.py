from scipy.stats import chi2_contingency
from sample_sizer import p_value, Variant
from expects import *
from unittest import TestCase

class TestSampleSizeCalculator(TestCase):
    def test_variant_can_check_unconverted(self):
        variant = Variant(conversions=1000,population=2000)
        expect(variant.unconverted()).to(equal(1000))

    def test_variant_can_check_conversion_rate(self):
        variant = Variant(conversions=1000,population=2000)
        expect(variant.conversion_rate()).to(equal(0.5))

    def test_can_check_p_value(self):
        control = Variant(conversions=11,population=43)
        variant = Variant(conversions=12,population=24)

        p_val = p_value(
            control = control,
            variant = variant
        )
        exp_t_statistic, exp_pval , dof , exp = chi2_contingency([[12,12],[11,32]], correction=False)
        expect(round(p_val, 6)).to(equal(round(exp_pval,6)))