from dataclasses import dataclass
from scipy.stats import chi2, chi2_contingency

@dataclass
class Variant:
    """Class to handle simple calculations for a variant"""
    conversions: int
    population: int

    def unconverted(self) --> int:
        return self.population - self.conversions

    def conversion_rate(self) --> float:
        return self.conversions/self.population

def p_value(control:Variant, variant:Variant) -> float:
    """
    Assumption: Conversions in control and variant follow the same binomial distribution B(1,p)

              Converted  Unconverted      Total
      Control   100           900         1000
      Variant   250           1750        2000
      Total     350           2650        3000

    """
    total_conversions = control.conversions + variant.conversions
    total_unconverted = control.unconverted() + variant.unconverted()
    total_population = control.population + variant.population

    exp_control_conversions = control.population*(total_conversions/total_population)
    exp_variant_conversions = variant.population*(total_conversions/total_population)
    exp_control_unconversions = control.population*(total_unconverted/total_population)
    exp_variant_unconversions = variant.population*(total_unconverted/total_population)

    test_statistic = (
        (abs(control.conversions - exp_control_conversions))**2/exp_control_conversions + 
        (abs(variant.conversions - exp_variant_conversions))**2/exp_variant_conversions +
        (abs(control.unconverted() - exp_control_unconversions))**2/exp_control_unconversions +
        (abs(variant.unconverted() - exp_variant_unconversions))**2/exp_variant_unconversions
    )

    degrees_of_freedom = 1
    p_value = 1 - chi2.cdf(test_statistic, degrees_of_freedom)
    return p_value
    