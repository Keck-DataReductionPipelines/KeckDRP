

class KcwiImageTypes(RecipesEngine.FrameType):
    def __init__(self):
        super(KcwiImageTypes, self).__init__()


KcwiImageTypes = KcwiImageTypes()


# bias: IMGTYPE = "BIAS" on extension 1

bias = [{'extension': 1, 'keyword': 'IMGTYPE', 'value': 'BIAS'}]

KcwiImageTypes.add_type(bias,'BIAS')
