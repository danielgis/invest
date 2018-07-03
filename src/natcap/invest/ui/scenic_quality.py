# coding=UTF-8

from natcap.invest.ui import model, inputs
from natcap.invest.scenic_quality import scenic_quality


class ScenicQuality(model.InVESTModel):
    def __init__(self):
        model.InVESTModel.__init__(
            self,
            label=u'Scenic Quality',
            target=scenic_quality.execute,
            validator=scenic_quality.validate,
            localdoc=u'../documentation/scenic_quality.html')

        self.general_tab = inputs.Container(
            interactive=True,
            label=u'General')
        self.add_input(self.general_tab)
        self.aoi_path = inputs.File(
            args_key=u'aoi_path',
            helptext=(
                u"An OGR-supported vector file.  This AOI instructs "
                u"the model where to clip the input data and the extent "
                u"of analysis.  Users will create a polygon feature "
                u"layer that defines their area of interest.  The AOI "
                u"must intersect the Digital Elevation Model (DEM)."),
            label=u'Area of Interest (Vector) (Required)',
            validator=self.validator)
        self.general_tab.add_input(self.aoi_path)
        self.structure_path = inputs.File(
            args_key=u'structure_path',
            helptext=(
                u"An OGR-supported vector file.  The user must specify "
                u"a point feature layer that indicates locations of "
                u"objects that contribute to negative scenic quality, "
                u"such as aquaculture netpens or wave energy "
                u"facilities.  In order for the viewshed analysis to "
                u"run correctly, the projection of this input must be "
                u"consistent with the project of the DEM input."),
            label=u'Features Impacting Scenic Quality (Vector) (Required)',
            validator=self.validator)
        self.general_tab.add_input(self.structure_path)
        self.dem_path = inputs.File(
            args_key=u'dem_path',
            helptext=(
                u"A GDAL-supported raster file.  An elevation raster "
                u"layer is required to conduct viewshed analysis. "
                u"Elevation data allows the model to determine areas "
                u"within the AOI's land-seascape where point features "
                u"contributing to negative scenic quality are visible."),
            label=u'Digital Elevation Model (Raster) (Required)',
            validator=self.validator)
        self.general_tab.add_input(self.dem_path)
        self.refraction = inputs.Text(
            args_key=u'refraction',
            helptext=(
                u"The earth curvature correction option corrects for "
                u"the curvature of the earth and refraction of visible "
                u"light in air.  Changes in air density curve the light "
                u"downward causing an observer to see further and the "
                u"earth to appear less curved.  While the magnitude of "
                u"this effect varies with atmospheric conditions, a "
                u"standard rule of thumb is that refraction of visible "
                u"light reduces the apparent curvature of the earth by "
                u"one-seventh.  By default, this model corrects for the "
                u"curvature of the earth and sets the refractivity "
                u"coefficient to 0.13."),
            label=u'Refractivity Coefficient (Required)',
            validator=self.validator)
        self.general_tab.add_input(self.refraction)
        self.valuation_tab = inputs.Container(
            interactive=True,
            label=u'Valuation')
        self.add_input(self.valuation_tab)
        self.valuation_function = inputs.Dropdown(
            args_key=u'valuation_function',
            helptext=(
                u"This field indicates the functional form f(x) the "
                u"model will use to value the visual impact for each "
                u"viewpoint.  For distances less than 1 km (x<1), the "
                u"model uses a linear form g(x) where the line passes "
                u"through f(1) (i.e.  g(1) == f(1)) and extends to zero "
                u"with the same slope as f(1) (i.e.  g'(x) == f'(1))."),
            label=u'Valuation Function',
            options=[u'polynomial: a + bx + cx^2 + dx^3',
                     u'logarithmic: a + b log(x)',
                     u'exponential: a * e^(-bx)'])
        self.valuation_tab.add_input(self.valuation_function)
        self.a_coefficient = inputs.Text(
            args_key=u'a_coef',
            helptext=(
                u"First coefficient used either by the polynomial or "
                u"by the logarithmic valuation function."),
            label=u"'a' Coefficient (polynomial/logarithmic/exponential)",
            validator=self.validator)
        self.valuation_tab.add_input(self.a_coefficient)
        self.b_coefficient = inputs.Text(
            args_key=u'b_coef',
            helptext=(
                u"Second coefficient used either by the polynomial or "
                u"by the logarithmic valuation function."),
            label=u"'b' Coefficient (polynomial/logarithmic/exponential)",
            validator=self.validator)
        self.valuation_tab.add_input(self.b_coefficient)
        self.c_coefficient = inputs.Text(
            args_key=u'c_coef',
            helptext=u"Third coefficient for the polynomial's quadratic term.",
            label=u"'c' Coefficient (polynomial only)",
            validator=self.validator)
        self.valuation_tab.add_input(self.c_coefficient)
        self.d_coefficient = inputs.Text(
            args_key=u'd_coef',
            helptext=u"Fourth coefficient for the polynomial's cubic exponent.",
            label=u"'d' Coefficient (polynomial only)",
            validator=self.validator)
        self.valuation_tab.add_input(self.d_coefficient)
        self.max_valuation_radius = inputs.Text(
            args_key=u'max_valuation_radius',
            helptext=(
                u"Radius beyond which the valuation is set to zero. "
                u"The valuation function 'f' cannot be negative at the "
                u"radius 'r' (f(r)>=0)."),
            label=u'Maximum Valuation Radius (meters) (Required)',
            validator=self.validator)
        self.valuation_tab.add_input(self.max_valuation_radius)

    def assemble_args(self):
        args = {
            self.workspace.args_key: self.workspace.value(),
            self.suffix.args_key: self.suffix.value(),
            self.aoi_path.args_key: self.aoi_path.value(),
            self.structure_path.args_key: self.structure_path.value(),
            self.dem_path.args_key: self.dem_path.value(),
            self.refraction.args_key: self.refraction.value(),
            self.valuation_function.args_key: self.valuation_function.value(),
            self.a_coefficient.args_key: self.a_coefficient.value(),
            self.b_coefficient.args_key: self.b_coefficient.value(),
            self.c_coefficient.args_key: self.c_coefficient.value(),
            self.d_coefficient.args_key: self.d_coefficient.value(),
            self.max_valuation_radius.args_key:
                self.max_valuation_radius.value(),
        }

        return args
