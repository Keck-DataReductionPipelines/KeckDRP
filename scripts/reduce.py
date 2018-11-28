#!/usr/bin/env python
import KeckDRP
import argparse
import importlib
import os
import glob
import time
import sys
from KeckDRP import conf
from KeckDRP import Instruments
from astropy import log


log.setLevel('INFO')

parser = argparse.ArgumentParser(description="""Perform a reduction.""",
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--recipe', type=str, help='reduction recipe')
parser.add_argument('--loop', action='store_true', help='Use infinite loop')
parser.add_argument('--imlist',
                    help='File containing the frames to be reduced')
parser.add_argument('frame', nargs='?', type=str, help='input image file')

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



    # load the frame and instantiate the object
    if os.path.isfile(image):
        frame = KeckDRP.KcwiCCD.read(image, unit='adu')
    else:
        log.error("The specified file (%s) does not exist" % image)
        sys.exit(1)

    inst = find_instrument(frame)

    if inst=='KCWI':
        Instrument = Instruments.KCWI()

    frame_type = Instrument.get_image_type(frame)

    recipe = Instrument.get_recipe(frame_type)

#    try:
#        mymodule = importlib.import_module(f'KeckDRP.{inst}.recipes.{recipe}')
    mymodule = importlib.import_module("KeckDRP.%s.recipes.%s" % (inst, recipe))
    myrecipe = getattr(mymodule, recipe)
#    except:
#        log.warn(f'\n--- Recipe {recipe} does not exist')
#        log.warn("\n--- Recipe %s does not exist" % (recipe))
#        return

    log.info("\n---  Reducing frame %s with recipe: %s ---" % (image, myrecipe.__name__))
    p = Instrument.get_primitives_class()
    myrecipe(p, frame)


def check_redux_dir():
    if not os.path.isdir(conf.REDUXDIR):
        os.makedirs(conf.REDUXDIR)
        log.info("Output directory created: %s" % conf.REDUXDIR)


def find_instrument(frame):
    if 'KCWI' in frame.header['INSTRUME']:
        return 'KCWI'


if __name__ == '__main__':

    args = parser.parse_args()

    check_redux_dir()

    if args.frame:
        log.info("reducing image %s" % args.frame)
        main_loop(frame=args.frame, recipe=args.recipe)
    elif args.loop:
        log.info("reducing in a loop")
        main_loop(loop=args.loop, recipe=args.recipe)
    elif args.imlist:
        log.info("reducing images in list in %s" % args.imlist)
        main_loop(imlist=args.imlist, recipe=args.recipe)
    else:
        log.info("Must supply an image or a list, or loop")
