from .GenomeTrack import GenomeTrack
from ..utilities import opener, count_lines
import numpy as np
from intervaltree import IntervalTree, Interval
from tqdm import tqdm

from .. readGwas import ReadGwas

DEFAULT_GWAS_COLOR = '#ff7f00'


class GwasTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.gwas', '.linear', '.logistic', '.assoc', '.qassoc']  # this is used by make_tracks_file to guess the type of track based on file name
    TRACK_TYPE = 'gwas'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# File containing the data. We expect an IGV .gwas format file with the columns: CHR, BP, SNP and P.
# Optionally, extra annotation columns can be added.
file =
# Indicate if your file has a header:
file_has_header = false
# Each SNP will be plotted as a 'o' and you can control color/size etc...
# Inside color
#color = red
# Border color
#border_color = black
# Line width
#line_width = 0.5
# Size
#marker_size = 45
# Optional. If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {'orientation': None,
                           'color': DEFAULT_GWAS_COLOR,
                           'border_color': 'black',
                           'line_width': 0.5,
                           'marker_size': 45,
                           'file_has_header': False}

    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = ['file_has_header']
    STRING_PROPERTIES = ['title', 'file_type', 'file', 'color', 'border_color']
    FLOAT_PROPERTIES = {'height': [0, np.inf],
                        'marker_size': [0, np.inf],
                        'line_width': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def __init__(self, *args, **kwarg):
        super(GwasTrack, self).__init__(*args, **kwarg)
        self.interval_tree = self.process_gwas(self.properties['region'])

    def process_gwas(self, plot_regions=None):
        """Read the gwas file and store values in a IntervalTree

        :param list plot_regions: list of plotted regions (like [(chrom1, start1, end1), (chrom2, start2, end2)]), defaults to None
        :return None
        """

        total_length = count_lines(opener(self.properties['file']),
                                   asBed=True)
        gwas_file_h = ReadGwas(opener(self.properties['file']),
                               has_header=self.properties['file_has_header'])

        valid_intervals = 0
        interval_tree = {}

        if plot_regions is not None:
            chroms_to_plot = set([v[0] for v in plot_regions])
        else:
            chroms_to_plot = None

        for record in tqdm(gwas_file_h, total=total_length):

            if plot_regions is not None and record.chromosome not in chroms_to_plot:
                continue
            
            if record.chromosome not in interval_tree:
                interval_tree[record.chromosome] = IntervalTree()

            interval_tree[record.chromosome].add(Interval(record.position,
                                                          record.position + 1, record))
            valid_intervals += 1

        try:
            gwas_file_h.file_handle.close()
        except AttributeError:
            pass

        if valid_intervals == 0:
            self.log.warning("No valid intervals were found in file "
                             f"{self.properties['file']} for regions"
                             f"{plot_regions}.\n")

        return interval_tree

    def plot(self, ax, chrom_region, start_region, end_region):
        """
        Plot a scatter plot for the GWAS data.
        The p-values are transformed as -log10(pvalue), so the y-axis will show the exponents of the p-values.

        :param ax: matplotlib axis
        :param chrom_region: chromosome name
        :param start_region: start position of the region
        :param end_region: end position of the region
        :return: None
        """
        if chrom_region not in self.interval_tree.keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.interval_tree.keys():
                self.log.warning("*Warning*\nNo interval was found when "
                                 "overlapping with both "
                                 f"{chrom_region_before}:{start_region}-{end_region}"
                                 f" and {chrom_region}:{start_region}-{end_region}"
                                 " inside the gwas file. "
                                 "This will generate an empty track!!\n")
                return

        gwas_overlap = \
            self.interval_tree[chrom_region][start_region:end_region]

        # Fill in the position and pvalues lists with data from the GWAS file
        position = [region.begin for region in gwas_overlap]
        # Notice the -log10 transformation
        y_values = [-np.log10(region.data.pvalue) if region.data.pvalue > 0 else 0
                    for region in gwas_overlap]

        # Plot the scatterplot
        ax.scatter(position, y_values,
                   s=self.properties['marker_size'],
                   color=self.properties['color'], marker='o',
                   edgecolors=self.properties['border_color'],
                   linewidths=self.properties['line_width'])
