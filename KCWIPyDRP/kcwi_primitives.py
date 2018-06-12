from KCWIPyDRP.ccd import ccd
from KCWIPyDRP.imgmath import imgmath

class kcwi(ccd,imgmath):
    def __init__(self):
        super(kcwi,self).__init__()
        print("init")
    def output_master(self,master_type="BIAS"):
        print("output_master %s" % master_type)
    def subtract_scattered_light(self):
        print("subtract_scattered_light")

