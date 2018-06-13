from KCWIPyDRP import kcwi_primitives
import argparse
import importlib
import os
from KCWIPyDRP.kcwi_objects import KcwiCCD
from astropy import log
global log
log.setLevel('INFO')


def go(iimg, rcp):
    if rcp is None:
        log.info("Checking %s header for recipe" % iimg)
        quit()

    mymodule = importlib.import_module("KCWIPyDRP.recipes."+str(rcp))
    recipe = getattr(mymodule, rcp)

    # load the frame and instantiate the object
    if os.path.isfile(iimg):
        global frame
        frame = KcwiCCD.read(iimg, unit='adu')
    else:
        log.info("The specified file name does not exist")
        sys.exit(1)

    log.info("executing recipe: %s" % rcp)
    p = kcwi_primitives.KcwiPrimitives()
    recipe(p, frame)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """Perform a reduction.
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-recipe', type=str, help='reduction recipe')
    parser.add_argument('inimg', type=str, help='input image file')

    args = parser.parse_args()

    log.info("input  image: %s" % args.inimg)
    log.info("input recipe: %s" % args.recipe)

    go(args.inimg, args.recipe)

