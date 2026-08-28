# -*- coding: utf-8 -*-
import os.path
from tempfile import NamedTemporaryFile

import matplotlib as mpl
from get_matplotlib_CI_version import get_CI_mpl_version
from matplotlib.testing.compare import compare_images

import pygenometracks.plotTracks

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks = """
[fixed line width]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = fixed link_line_width to 0.5 link_fontsize = 10 orientation = inverted
height = 3
bw_color = darkblue
max_value = 110
nans_to_zeros = true
show_data_range = true
link_color = darkblue
link_fontsize = 10
link_line_width = 0.5
orientation = inverted
file_type = sashimiBigWig

[no label]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = link_labels = false; alpha = 0.1; link_scale_height = 10
height = 3
bw_color = purple
alpha = 0.1
max_value = 110
show_data_range = true
link_color = purple
link_labels = false
link_scale_height = 10
file_type = sashimiBigWig

[spacer]
height = 0.1

[scale line width]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = link_scale_line_width = 10; link_fontsize = 4; link_alpha = 0.5
height = 3
bw_color = orange
max_value = 110
show_data_range = true
link_color = orange
link_fontsize = 4
link_alpha = 0.5
link_scale_line_width = 10
file_type = sashimiBigWig

[log transformed]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = transform = log; log_pseudocount = 1
height = 3
bw_color = black
max_value = 110
nans_to_zeros = true
show_data_range = true
link_color = black
transform = log
log_pseudocount = 1
file_type = sashimiBigWig

[log transformed min_value]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = transform = log; log_pseudocount = 1; min_value = 50
height = 3
bw_color = black
min_value = 50
max_value = 110
nans_to_zeros = true
show_data_range = true
link_color = black
transform = log
log_pseudocount = 1
file_type = sashimiBigWig

[log transformed min_value original y_axis_values]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = transform = log; log_pseudocount = 1; min_value = 10; y_axis_values = original
height = 3
bw_color = grey
min_value = 10
max_value = 110
nans_to_zeros = true
show_data_range = true
link_color = grey
transform = log
log_pseudocount = 1
y_axis_values = original
file_type = sashimiBigWig


[min_value grid]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = min_value = 10; grid = true; link_color = viridis; link_line_width = 2
height = 3
bw_color = cyan
min_value = 10
max_value = 110
nans_to_zeros = true
show_data_range = true
link_line_width = 2
link_color = viridis
grid = true
file_type = sashimiBigWig

[colormap min_value max_value]
file = avg_chr2-231091223_231109786_231113600_0.bw
link_file = chr2-231091223_231109786_231113600_0.sashimi
title = link_color = viridis; link_line_width = 2; link_min_max_values to 0 1
height = 3
bw_color = pink
max_value = 110
show_data_range = no
link_line_width = 2
link_color = viridis
link_min_value = 0
link_max_value = 1
link_labels = false
link_scale_height = 10
file_type = sashimiBigWig
"""

with open(os.path.join(ROOT, "sashimi_tracks.ini"), 'w') as fh:
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance
default_mpl_version = get_CI_mpl_version()


def test_sashimi_main():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 29
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "sashimi_tracks.ini")
    region = "chr2:231107879-231115507"
    expected_file = os.path.join(ROOT, 'master_sashimi.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_sashimi_X():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 28
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "sashimi_tracks.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_sashimi_X.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)
