from astropy.io import fits
import glob
import shutil


def get_keyword(hdu_in, keyword):
    keyword = hdu_in[0].header[keyword]
    if type(keyword) is str:
        keyword = keyword.strip()
    return keyword


files = glob.glob("kb*.fits")
new_files = sorted(files)

group_id = ""
previous_ccdcfg = ""
previous_frameno = -1
previous_groupid = ""
previous_imtype = ""
previous_stateid = ""
previous_ttime = -1.

for file in new_files:
    hdu = fits.open(file, mode='update')
    shutil.copy(file, "%s.bak" % file)
    camera = get_keyword(hdu, 'CAMERA')
    if 'CCDCFG' not in hdu[0].header:
        ccdcfg = get_keyword(hdu, 'CCDSUM').replace(" ", "")
        ccdcfg += "%1d" % get_keyword(hdu, 'CCDMODE')
        ccdcfg += "%02d" % get_keyword(hdu, 'GAINMUL')
        ccdcfg += "%02d" % get_keyword(hdu, 'AMPMNUM')
        hdu[0].header['CCDCFG'] = ccdcfg
    if file == new_files[0]:
        # first file
        current_imtype = get_keyword(hdu, 'IMTYPE')
        current_frameno = get_keyword(hdu, 'FRAMENO')
        current_stateid = get_keyword(hdu, 'STATEID')
        current_ttime = get_keyword(hdu, 'TTIME')
        current_ccdcfg = get_keyword(hdu, 'CCDCFG')
        group_id = "%s-%d" % (get_keyword(hdu, 'DATE-OBS'),
                              get_keyword(hdu, 'FRAMENO'))
        print("File: %s, ImType: %s, GroupId: %s" % (file, current_imtype,
                                                     group_id))
        previous_imtype = current_imtype
        previous_frameno = current_frameno
        previous_stateid = current_stateid
        previous_ttime = current_ttime
        previous_groupid = group_id
        previous_ccdcfg = current_ccdcfg

    else:
        # other files
        current_imtype = get_keyword(hdu, 'IMTYPE')
        current_frameno = get_keyword(hdu, 'FRAMENO')
        current_ccdcfg = get_keyword(hdu, 'CCDCFG')
        current_stateid = get_keyword(hdu, 'STATEID')
        current_ttime = get_keyword(hdu, 'TTIME')

        # bias
        if current_imtype == 'BIAS' and previous_imtype == 'BIAS':
            if current_frameno == previous_frameno + 1 and \
                    current_ccdcfg == previous_ccdcfg:
                group_id = previous_groupid

        # flat
        elif current_imtype == 'FLATLAMP' and previous_imtype == "FLATLAMP":
            if current_frameno == previous_frameno + 1 and \
                    current_stateid == previous_stateid:
                group_id = previous_groupid

        # domeflat
        elif current_imtype == 'DOMEFLAT' and previous_imtype == "DOMEFLAT":
            if current_frameno == previous_frameno + 1 and \
                    current_stateid == previous_stateid:
                group_id = previous_groupid

        # twilight flight
        elif current_imtype == 'TWIFLAT' and previous_imtype == "TWIFLAT":
            if current_frameno == previous_frameno + 1 and \
                    current_stateid == previous_stateid:
                group_id = previous_groupid

        # darks
        elif current_imtype == 'DARK' and previous_imtype == "DARK":
            if current_frameno == previous_frameno + 1 and \
                    current_stateid == previous_stateid and \
                    current_ttime == previous_ttime:
                group_id = previous_groupid

        else:
            # need a new group id
            group_id = "%s-%d" % (get_keyword(hdu, 'DATE-OBS'),
                                  get_keyword(hdu, 'FRAMENO'))

        print("File: %s, ImType: %s, GroupId: %s" % (file, current_imtype,
                                                     group_id))

        hdu[0].header['GROUPID'] = group_id
        hdu.flush()
        hdu.close()
        previous_imtype = current_imtype
        previous_frameno = current_frameno
        previous_ccdcfg = current_ccdcfg
        previous_stateid = current_stateid
        previous_groupid = group_id
        previous_ttime = current_ttime
