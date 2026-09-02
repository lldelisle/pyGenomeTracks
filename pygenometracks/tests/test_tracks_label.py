# -*- coding: utf-8 -*-
import os.path
from tempfile import NamedTemporaryFile, TemporaryDirectory

import matplotlib as mpl
from get_matplotlib_CI_version import get_CI_mpl_version
from matplotlib.testing.compare import compare_images

import pygenometracks.plotTracks

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks = """
[x-axis]
title = A very long title that does not fit into a single line whatever is the width of the label
"""
with open(os.path.join(ROOT, "title.ini"), 'w') as fh:
    fh.write(browser_tracks)


tolerance = 13  # default matplotlib pixed difference tolerance
default_mpl_version = get_CI_mpl_version()


def test_regular_width_label():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 51
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.2.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 53
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_ral():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 50
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_ral.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           "--trackLabelHAlign right "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_cal():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 50
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_cal.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           "--trackLabelHAlign center "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_cal_dpi250():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 54
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_cal_d250.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 250 "\
           "--trackLabelHAlign center "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_big_font():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 92
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_fs20.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           "--fontSize 20 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fixed_height():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 24
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_force_height.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--height 10 --title force_height "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_non_existing_dir():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 51
    else:
        my_tolerance = tolerance

    outdir = TemporaryDirectory()
    output_file = os.path.join(outdir.name, "pGT_test", "test.png")
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.2.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {output_file}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         output_file, my_tolerance)
    assert res is None, res

    outdir.cleanup()

def test_bed_with_multiple_regions_with_identical_names():

    ini_file = os.path.join(ROOT, "title.ini")

    outdir = TemporaryDirectory()
    bed_file = os.path.join(outdir.name, 'test.bed')
    with open(bed_file, 'w') as f:
        f.write('X\t0\t10\tname1\n')
        f.write('X\t0\t20\tname1\n')
        f.write('X\t0\t50\tname2\n')
    output_file = os.path.join(outdir.name, "test.png")
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {output_file}".split()
    pygenometracks.plotTracks.main(args)
    all_files = os.listdir(outdir.name)
    assert len(all_files) == 4
    assert set(all_files) == {'test.bed', 'test_name1.png', 'test_name2.png', 'test_X-0-20.png'}
    