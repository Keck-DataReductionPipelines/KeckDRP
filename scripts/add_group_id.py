from astropy.io import fits
import glob
import shutil

def get_keyword(hdu, keyword):
    keyword =  hdu[0].header[keyword]
    if type(keyword) is str:
        keyword = keyword.strip()
    return keyword
def assign_group_id(file, group_id, imtype):
    print("File: %s, ImType: %s, GroupId: %s" % (file, imtype, group_id))


files = glob.glob("kb*.fits")
new_files = sorted(files)

for file in new_files:
    hdu = fits.open(file, mode='update')
    shutil.copy(file, "%s.bak" % file)
    camera = get_keyword(hdu, 'CAMERA')
    if file == new_files[0]:
        # print("first file")
        current_imtype = get_keyword(hdu, 'IMTYPE')
        current_frameno = get_keyword(hdu, 'FRAMENO')
        current_stateid = get_keyword(hdu, 'STATEID')
        current_ttime = get_keyword(hdu, 'TTIME')
        current_ccdcfg = get_keyword(hdu, 'CCDCFG')
        group_id = "%s-%d" % (get_keyword(hdu, 'DATE-OBS'), get_keyword(hdu, 'FRAMENO'))
        assign_group_id(file, group_id, current_imtype)
        previous_imtype = current_imtype
        previous_frameno = current_frameno
        previous_stateid = current_stateid
        previous_ttime = current_ttime
        previous_groupid = group_id
        previous_ccdcfg = current_ccdcfg


    else:
        # print("other files")
        current_imtype = get_keyword(hdu, 'IMTYPE')
        current_frameno = get_keyword(hdu, 'FRAMENO')
        current_ccdcfg = get_keyword(hdu, 'CCDCFG')
        current_stateid = get_keyword(hdu, 'STATEID')
        current_ttime = get_keyword(hdu, 'TTIME')
        # print("Working on file: %s (frameno: %d, imtype: %s)" % (file,current_frameno, current_imtype))

        # bias

        if current_imtype == 'BIAS' and previous_imtype == 'BIAS':
            if current_frameno == previous_frameno + 1 and current_ccdcfg == previous_ccdcfg:
                group_id = previous_groupid

        # flat

        elif current_imtype == 'FLATLAMP' and previous_imtype == "FLATLAMP":
            if current_frameno == previous_frameno + 1 and current_stateid == previous_stateid:
                group_id = previous_groupid

        # domeflat

        elif current_imtype == 'DOMEFLAT' and previous_imtype == "DOMEFLAT":
            if current_frameno == previous_frameno + 1 and current_stateid == previous_stateid:
                group_id = previous_groupid

        # twilight flight

        elif current_imtype == 'TWIFLAT' and previous_imtype == "TWIFLAT":
            if current_frameno == previous_frameno + 1 and current_stateid == previous_stateid:
                group_id = previous_groupid

        # darks

        elif current_imtype == 'DARK' and previous_imtype == "DARK":
            if current_frameno == previous_frameno + 1 and current_stateid == previous_stateid and current_ttime == prevoiius_ttime:
                group_id = previous_groupid

        else:
            # need a new group id
            group_id = "%s-%d" % (get_keyword(hdu, 'DATE-OBS'), get_keyword(hdu, 'FRAMENO'))

        assign_group_id(file, group_id, current_imtype)
        if "BLUE" in camera:
            keyword = "GROUP-B"
        elif "RED" in camera:
            keyword = "GROUP-R"
        else:
            print("No camera information, stopping")
            continue
        hdr = hdu[0].header
        hdr[keyword] = group_id
        hdu.flush()
        hdu.close()
        previous_imtype = current_imtype
        previous_frameno = current_frameno
        previous_ccdcfg = current_ccdcfg
        previous_stateid = current_stateid
        previous_groupid = group_id
        previous_ttime = current_ttime


