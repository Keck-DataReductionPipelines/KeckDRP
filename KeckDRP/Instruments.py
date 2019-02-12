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
                actual_value = frame.header[keyword].strip()
                if actual_value != expected_value:
                    current_type_is_true = False
            if current_type_is_true is False:
                continue
            else:
                return type
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
            bias=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'BIAS'}],
            flatlamp=[{'extension': 1, 'keyword': 'IMTYPE',
                       'value': 'FLATLAMP'}],
            domeflat=[{'extension': 1, 'keyword': 'IMTYPE',
                       'value': 'DOMEFLAT'}],
            contbars=[{'extension': 1, 'keyword': 'IMTYPE',
                       'value': 'CONTBARS'}],
            arclamp=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'ARCLAMP'}],
            object=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'OBJECT'}])
        self.recipes = dict(
            focus=None,
            test=None,
            bias='proccess_biases',
            flatlamp='handle_flats',
            domeflat='make_master_dome_flat',
            contbars='handle_contbars',
            arclamp='handle_arclamp',
            object='make_science')

    def get_primitives_class(self):
        from KeckDRP.KCWI import kcwi_primitives
        return kcwi_primitives.KcwiPrimitives()
