"""Keck Instrument Classes"""


class Instrument:
    """Define generic instrument class properties and methods

    Here we define the generic init function, a way to derive the
    image type, how to get the reduction recipe based on the image type
    and the method for accessing primitives.
    """

    def __init__(self):
        self.image_types = {}
        self.recipes = {}

    def get_image_type(self, frame):
        # loop through the different image types
        for image_type in self.image_types.keys():
            current_type_is_true = True
            for keyword_check in self.image_types[image_type]:
                extension = keyword_check['extension']
                keyword = keyword_check['keyword']
                expected_value = keyword_check['value']
                actual_value = frame.header[keyword]
                if type(actual_value) is str:
                    actual_value = actual_value.strip()
                if actual_value != expected_value:
                    current_type_is_true = False
            if current_type_is_true is False:
                continue
            else:
                return image_type
        return "UNKNOWN"

    def get_recipe(self, image_type):
        try:
            return self.recipes[image_type]
        except:
            return "None"

    def get_primitives_class(self):
        class generic_primitives():
            def __init__(self):
                pass


class KCWI(Instrument):
    """Define the KCWI Instrument Class

    Here we provide a dictionary that defines image types and how those
    types are mapped to reduction recipes.  We also override the generic
    primitive class with KCWI-specific primitives.
    """

    def __init__(self):
        super(KCWI, self).__init__()
        self.image_types = dict(
            focus=[{'extension': 1, 'keyword': 'OBJECT', 'value': 'focus'}],
            test=[{'extension': 1, 'keyword': 'OBJECT', 'value': 'TEST'}],
            bias=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'BIAS'},
                  {'extension': 1, 'keyword': 'TELAPSE', 'value': 0.0}],
            flatlamp=[{'extension': 1, 'keyword': 'IMTYPE',
                       'value': 'FLATLAMP'}],
            domeflat=[{'extension': 1, 'keyword': 'IMTYPE',
                       'value': 'DOMEFLAT'}],
            twiflat=[{'extension': 1, 'keyword': 'IMTYPE',
                     'value': 'TWIFLAT'}],
            contbars=[{'extension': 1, 'keyword': 'IMTYPE',
                       'value': 'CONTBARS'}],
            arclamp=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'ARCLAMP'}],
            object=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'OBJECT'}])
        self.recipes = dict(
            focus=None,
            test=None,
            bias='process_biases',
            dark='process_darks',
            flatlamp='process_internal_flats',
            domeflat='process_dome_flats',
            contbars='process_contbars',
            arclamp='process_arclamps',
            object='make_science')

    def get_primitives_class(self):
        from KeckDRP.KCWI import kcwi_primitives
        return kcwi_primitives.KcwiPrimitives()
