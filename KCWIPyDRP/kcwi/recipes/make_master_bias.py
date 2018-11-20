;
; Copyright (c) 2017, Chris Curtin, no rights reserved.
;+
; NAME:
;	curtin_make_master_bias.py
;
; PURPOSE:
;	This recipe should eventually consist of the primitives necessary
;   to receive as input a set of biases taken with Keck KCWI and output
;   a master bias of the specifications of the KCWI DRP to be used for
;   bias subtraction.
;
; CALLING SEQUENCE:
;	python make_master_bias.py p frame
;
; OPTIONAL INPUTS:
;   ____________________________ (unknown at this time)
;
; KEYWORDS:
;   ____________________________ (unknown at this time)
;
; OUTPUTS:
;	____________________________ (the naming convention of the master bias file)
;
;;;;;;;;;;;;;;;;;;;;;;;; SIDE EFFECTS:
::::::::::::::::::::::::	Outputs processed files in output directory specified by the
;;;;;;;;;;;;;;;;;;;;;;;;	KCWI_PPAR struct read in from Pparfname.
;   ____________________________ (unknown at this time, but follow the convention above)
;
;;;;;;;;;;;;;;;;;;;;;;;; PROCEDURE:
;;;;;;;;;;;;;;;;;;;;;;;;	Reads Pparfname to derive input/output directories and reads the
;;;;;;;;;;;;;;;;;;;;;;;;	corresponding '*.proc' file in output directory to derive the list
;;;;;;;;;;;;;;;;;;;;;;;;	of input bias files.  Each input
;;;;;;;;;;;;;;;;;;;;;;;;	file is read in.  The overscan region for each
;;;;;;;;;;;;;;;;;;;;;;;;	image is then analyzed and a row-by-row subtraction is performed
;;;;;;;;;;;;;;;;;;;;;;;;	to remove 1/f noise in the readout.  The images are trimmed and
;;;;;;;;;;;;;;;;;;;;;;;;	assembled into physical CCD-sized images.  Next, a gain correction
;;;;;;;;;;;;;;;;;;;;;;;;	for each amplifier is applied to convert each image into electrons.
;;;;;;;;;;;;;;;;;;;;;;;;   A master bias is generated from the dressed images
;;;;;;;;;;;;;;;;;;;;;;;;	which accounts for CCD read noise.
;
;;;;;;;;;;;;;;;;;;;;;;;; EXAMPLE:
;;;;;;;;;;;;;;;;;;;;;;;;	Perform stage1 reductions on the images in 'night1' directory and put
;;;;;;;;;;;;;;;;;;;;;;;;	results in 'night1/redux':
;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;	KCWI_PREP,'night1','night1/redux'
;;;;;;;;;;;;;;;;;;;;;;;;	KCWI_STAGE1,'night1/redux/kcwi.proc'
;
; MODIFICATION HISTORY:
;	Written by:	Chris Curtin (ccurtin@swin.edu.au) with the conventions
;   of Don Neill (neill@caltech.edu) and in collaboration with Luca Rizzi
;   (lrizzi@keck.hawaii.edu)
;	2018-NOV-15	Initial version
;-

def make_master_bias(p, frame):
    # attach frame data
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
    else:
        # update proc table
        p.update_proctab()
        p.write_proctab()
        p.log.info("bias counted")
        # how many biases do we have?
        combine_list = p.n_proctab(targtype='BIAS')
        p.log.info("number of biases = %d" % len(combine_list))
        # create master bias
        if len(combine_list) >= 7:
            p.img_combine(combine_list, keylog='BIASLIST')
            # output file and update proc table
            p.update_proctab(suffix='mbias', newtype='MBIAS')
            p.write_image(suffix='mbias')
            p.log.info("master bias produced")
        else:
            p.log.info("need 7 biases to produce master")
        p.write_proctab()