from KeckDRP import PrimitivesBASE
import KeckDRP
import ccdproc
import os

class DevelopmentPrimitives(PrimitivesBASE):

    def get_idl_counterpart(self, target_type=None):
        self.log.info("fetching idl counterpart")
        reference_tab = self.proctab[(self.proctab['OFNAME'] == self.frame.header['OFNAME'])]
        counterpart_tab = self.proctab[((self.proctab['TYPE'] == target_type)
                                        & (self.proctab['CID'] == reference_tab[0]['CID']))]
#        os.system('cp IDL_redux/'+counterpart_tab[0]['OFNAME']+'* .')
        return(counterpart_tab)

    def read_idl_copy(self, idl_reference_table=None, in_directory=None, suffix=None):
        self.log.info("accessing IDL redux")
        idl_image = idl_reference_table[0]['OFNAME']
        idl_number = idl_reference_table[0]['FRAMENO']

        if in_directory is None:
            prefix = 'IDL_redux/'
        else:
            prefix = in_directory

        if suffix is None:
            suffix = '.fits'
        else:
            suffix = '_' + suffix + '.fits'

        infile = os.path.join(prefix, idl_image.split('.')[0] + suffix)
        self.log.info("reading idl image: %s" % infile)
        idl_file = KeckDRP.KcwiCCD.read(infile, unit='adu')
        return(idl_file)


