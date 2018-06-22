#!/usr/bin/env python
from KCWIPyDRP.kcwi import kcwi_primitives
import argparse
import importlib
import os
import glob
import time
import sys
from KCWIPyDRP.kcwi.kcwi_objects import KcwiCCD
from astropy import log

log.setLevel('INFO')


def main_loop(frame=None, recipe=None, loop=None, imlist=None):
    # case 1: one one image is specified
    if frame:
        go(frame, recipe)
        return
    # case 2: infinite loop
    if loop:
        # step 1: build a list of files already in the directory
        input_files_before = glob.glob('k*.fits')
        for input_file in input_files_before:
            go(input_file, recipe)
        while True:
            # step 2: start infinite loop waiting for new files
            input_files_now = glob.glob('k*.fits')
            new_files = [f for f in input_files_now
                         if f not in input_files_before]
            for input_file in new_files:
                go(input_file, recipe)
            input_files_before = input_files_now
            time.sleep(10)
    # case 3: a list of files is provided
    if imlist:
        if os.path.isfile(imlist):
            with open(imlist) as file_list:
                for file in file_list:
                    go(file.rstrip('\n'), recipe)
        return


def go(image, rcp):
    if rcp is None:
        log.info("Checking %s header for recipe" % image)
        quit()

    mymodule = importlib.import_module("KCWIPyDRP.kcwi.recipes."+str(rcp))
    recipe = getattr(mymodule, rcp)

    # load the frame and instantiate the object
    if os.path.isfile(image):
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
    parser.add_argument('--imlist',
                        help='File containing the frames to be reduced')
    parser.add_argument('frame', nargs='?', type=str, help='input image file')

    args = parser.parse_args()

    if args.recipe is None:
        log.info("Must supply a recipe")
    else:
        if args.frame:
            log.info("applying recipe %s to image %s" % (args.recipe,
                                                         args.frame))
            main_loop(frame=args.frame, recipe=args.recipe)
        elif args.loop:
            log.info("applying recipe %s in a loop" % args.recipe)
            main_loop(recipe=args.recipe, loop=args.loop)
        elif args.imlist:
            log.info("applying recipe %s to img list in %s" % (args.recipe,
                                                               args.imlist))
            main_loop(recipe=args.recipe, imlist=args.imlist)
        else:
            log.info("Must supply an image or a list, or loop")

