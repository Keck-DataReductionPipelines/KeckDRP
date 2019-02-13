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


#os.system('rm -r redux')
#os.system('rm kcwi.proc')

log.setLevel('INFO')

parser = argparse.ArgumentParser(description="""Perform a reduction.""",
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--recipe', type=str, help='reduction recipe')
parser.add_argument('--loop', action='store_true', help='Use infinite loop')
parser.add_argument('--imlist',
                    help='File containing the frames to be reduced')
parser.add_argument('--imtype', type=str, help='reduce all frames of the specified image type')
parser.add_argument('--overwrite', action='store_true', help='Reprocess images, ignore proctab information')
parser.add_argument('frames', nargs='+', type=str, help='input image file')


def main_loop(frames=None, recipe=None, loop=None, imlist=None, imtype=None):
    # case 1: one one image is specified
    #if frames:
    #    frames_list = glob.glob(frames)
    for frame in frames:
        go(frame, recipe, imtype)
    return
    # case 2: infinite loop
    if loop:
        # step 1: build a list of files already in the directory
        input_files_before = glob.glob('*.fits')
        for input_file in input_files_before:
            go(input_file, recipe, imtype)
        while True:
            # step 2: start infinite loop waiting for new files
            input_files_now = glob.glob('*.fits')
            new_files = [f for f in input_files_now
                         if f not in input_files_before]
            for input_file in new_files:
                go(input_file, recipe, imtype)
            input_files_before = input_files_now
            time.sleep(10)
    # case 3: a list of files is provided
    if imlist:
        if os.path.isfile(imlist):
            with open(imlist) as file_list:
                for file in file_list:
                    go(file.rstrip('\n'), recipe, imtype)
        return


def go(image, rcp, imtype=None):

    # load the frame and instantiate the object
    if os.path.isfile(image):
        frame = KeckDRP.KcwiCCD.read(image, unit='adu')
        # handle missing CCDCFG
        if 'CCDCFG' not in frame.header:
            ccdcfg = frame.header['CCDSUM'].replace(" ", "")
            ccdcfg += "%1d" % frame.header['CCDMODE']
            ccdcfg += "%02d" % frame.header['GAINMUL']
            ccdcfg += "%02d" % frame.header['AMPMNUM']
            frame.header['CCDCFG'] = ccdcfg
    else:
        log.error("The specified file (%s) does not exist" % image)
        sys.exit(1)

    inst = find_instrument(frame)

    if inst == 'KCWI':
        Instrument = Instruments.KCWI()

    frame_type = Instrument.get_image_type(frame)
    if imtype:
        if imtype!=frame_type:
            log.info("Frame %s (%s) is not of imtype %s and will be skipped" % (image, frame_type, imtype))
            return

    recipe = Instrument.get_recipe(frame_type)
    log.info("Frame %s is of IMTYPE %s and uses the recipe %s" % (image, frame_type, recipe))

    if recipe is None:
        log.info("\n--- No reduction necessary ---\n")
        return

    try:
        mymodule = importlib.import_module("KeckDRP.%s.recipes.%s" % (inst, recipe))
        myrecipe = getattr(mymodule, recipe)
    except ImportError:
        log.warn("\n--- Recipe %s does not exist" % (recipe))
        return
    except RuntimeError:
        log.warn("\n--- Error looking for recipe %s" % (recipe))

    log.info("\n---  Reducing frame %s with recipe: %s ---" %
             (image, myrecipe.__name__))
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
    if args.overwrite:
        log.info("Reprocessing of file is enabled")
        conf.OVERWRITE = True

    if args.frames:
        log.info("reducing image(s) %s" % args.frames)
        main_loop(frames=args.frames, recipe=args.recipe, imtype=args.imtype)
    elif args.loop:
        log.info("reducing in a loop")
        main_loop(loop=args.loop, recipe=args.recipe, imtype=args.imtype)
    elif args.imlist:
        log.info("reducing images in list in %s" % args.imlist)
        main_loop(imlist=args.imlist, recipe=args.recipe, imtype=args.imtype)
    else:
        log.info("Must supply an image or a list, or loop")
