#!/usr/bin/env python
from KCWIPyDRP import kcwi_primitives
import argparse
import importlib
import os, glob, time, sys
from KCWIPyDRP.kcwi_objects import KcwiCCD
from astropy import log
global log
log.setLevel('INFO')


def main_loop(args):
    # case 1: one one image is specified
    if args.frame:
        go(args.frame, args.recipe)
        return
    # case 2: infinite loop
    if args.loop:
        # step 1: build a list of files already in the directory
        input_files_before = glob.glob('k*.fits')
        for input_file in input_files_before:
            go(input_file, args.recipe)
        while True:
            # step 2: start infinite loop waiting for new files
            input_files_now = glob.glob('k*.fits')
            new_files = [f for f in input_files_now if not f in input_files_before]
            for input_file in new_files:
                go(input_file, args.recipe)
            input_files_before = input_files_now
            time.sleep(10)
    # case 3: a list of files is provided
    if args.list:
        if os.path.isfile(args.list):
            with open(args.list) as file_list:
                for file in file_list:
                    go(file.rstrip('\n'), args.recipe)
        return


def go(image, rcp):
    if rcp is None:
        log.info("Checking %s header for recipe" % image)
        quit()

    mymodule = importlib.import_module("KCWIPyDRP.recipes."+str(rcp))
    recipe = getattr(mymodule, rcp)

    # load the frame and instantiate the object
    if os.path.isfile(image):
        global frame
        frame = KcwiCCD.read(image, unit='adu')
    else:
        log.info("The specified file (%s) does not exist" % (image))
        sys.exit(1)

    log.info("\n---  Reducing frame %s with recipe: %s ---" % (image, rcp))
    p = kcwi_primitives.KcwiPrimitives()
    recipe(p, frame)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """Perform a reduction.
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--recipe', type=str, help='reduction recipe')
    parser.add_argument('--loop', action='store_true', help='Use infinite loop')
    parser.add_argument('--list', help='File containing the frames to be reduced')
    parser.add_argument('frame', nargs='?',type=str, help='input image file')

    args = parser.parse_args()

    log.info("input  image: %s" % args.frame)
    log.info("input recipe: %s" % args.recipe)

    main_loop(args)

