import os.path
from tempfile import NamedTemporaryFile

import matplotlib as mpl
from get_matplotlib_CI_version import get_CI_mpl_version
from matplotlib.testing.compare import compare_images

import pygenometracks.plotTracks

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

tracks = """
[gwas]
file = gwas_1.gwas
height = 4
title = test_1 default values

[spacer]

[gwas_2]
file = gwas_2.gwas
file_has_header = True
height = 2
title = test_2 file_has_header = true color = #50E3C2 border_color = red line_width = 2 marker_size = 90 show_data_range = false
color = #50E3C2
border_color = red
line_width = 2
marker_size = 90
show_data_range = false

[spacer]

[gwas_2]
file = gwas_1.gwas
height = 4
title = test_1 default values min_value = 0 max_value = 15
min_value = 0
max_value = 15

[x-axis]
"""

with open(os.path.join(ROOT, "gwas.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance
default_mpl_version = get_CI_mpl_version()


def test_gwas_track():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 26
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='gwas_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gwas.ini")
    region = "X:3000000-3200000"
    expected_file = os.path.join(ROOT, 'master_gwas.png')
    args = f"--tracks {ini_file} --region {region} " \
           "--trackLabelFraction 0.2 --dpi 130 " \
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_gwas_track_chrX():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 15
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='gwas_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gwas.ini")
    region = "chrX:3000000-3200000"
    expected_file = os.path.join(ROOT, 'master_gwas.png')
    args = f"--tracks {ini_file} --region {region} " \
           "--trackLabelFraction 0.2 --dpi 130 " \
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance + 14)  # 14 corresponds to the 'chr' on the x axis
    assert res is None, res

    os.remove(outfile.name)


def test_gwas_track_chrY():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 26
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='gwas_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gwas.ini")
    region = "chrY:3000000-3200000"
    expected_file = os.path.join(ROOT, 'master_gwas_chrY.png')
    args = f"--tracks {ini_file} --region {region} " \
           "--trackLabelFraction 0.2 --dpi 130 " \
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)
