class Instrument():
    def __init__(self):
        self.image_types = {}
        self.recipes = {}

    def get_image_type(self, frame):
        # loop through the different image types
        for type in self.image_types.keys():
            current_type_is_true = True
            for keyword_check in self.image_types[type]:
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

    def get_recipe(self, type):
        try:
            return self.recipes[type]
        except:
            return "None"

    def get_primitives_class(self):
        class generic_primitives():
            def __init__(self):
                pass



class KCWI(Instrument):
    def __init__(self):
        super(KCWI, self).__init__()
        self.image_types = dict(bias=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'BIAS'}],
                           flatlamp=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'FLATLAMP'}],
                           domeflat=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'DOMEFLAT'}],
                           contbars=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'CONTBARS'}],
                           arclamp=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'ARCLAMP'}],
                           object=[{'extension': 1, 'keyword': 'IMTYPE', 'value': 'OBJECT'}])
        self.recipes = dict(bias='handle_biases',
                            flatlamp='make_master_flat',
                            domeflat='make_master_dome_flat',
                            contbars='reduce_contbars',
                            arclamp='reduce_arclamp',
                            object='reduce_science')
    def get_primitives_class(self):
        from KeckDRP.KCWI import kcwi_primitives
        return kcwi_primitives.KcwiPrimitives()


