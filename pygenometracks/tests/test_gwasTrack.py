import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

tracks = """
[gwas]
file = gwas_1.gwas
height = 4
title = test_1

[spacer]

[gwas_2]
file = gwas_2.gwas
file_has_header = True
height = 2
title = test_2
color = #50E3C2

[x-axis]
"""

with open(os.path.join(ROOT, "gwas.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_gwas_track():
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
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
