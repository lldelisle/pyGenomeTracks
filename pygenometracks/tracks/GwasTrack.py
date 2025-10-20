from .GenomeTrack import GenomeTrack
import numpy as np
from .. readGwas import ReadGwas

DEFAULT_GWAS_COLOR = '#ff7f00'


class GwasTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.gwas', '.linear', '.logistic', '.assoc', '.qassoc']  # this is used by make_tracks_file to guess the type of track based on file name
    TRACK_TYPE = 'gwas'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# Title of the track. Usually displayed to the right as a label
title =
# Height of the track    
height =
# File containing the data. We expect an IGV .gwas format file with the columns: CHR, BP, SNP and P. Optionally, extra
# annotation columns can be added.
file =
# Y label text
ylabel =
# Fontsize of the labels
fontsize =
# Color
color =
# Optional. If not given is guessed from the file ending.
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {'fontsize': 12,
                           'orientation': None,
                           'color': DEFAULT_GWAS_COLOR,
                           'border_color': 'black',
                           'labels': True,
                           'line_width': 0.5,
                           'max_labels': 60,
                           'max_value': 1,
                           'min_value': 0,
                           'fontstyle': 'normal',
                           'y_axis_max_val': None,
                           'id_fontsize': 12,
                           'marker_size': 45,
                           'file_has_header': False}

    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {}
    POSSIBLE_PROPERTIES = {}
    BOOLEAN_PROPERTIES = ['file_has_header']
    STRING_PROPERTIES = ['title', 'file_type', 'file', 'color']
    FLOAT_PROPERTIES = {'height': [0, np.inf], 'fontsize': [0, np.inf], 'id_fontsize': [0, np.inf],
                        'marker_size': [0, np.inf], 'y_axis_max_val': [0, np.inf]}
    INTEGER_PROPERTIES = {}

    def plot(self, ax, chrom, region_start, region_end):
        """
        Plot a scatter plot for the GWAS data.
        The p-values are transformed as -log10(pvalue), so the y-axis will show the exponents of the p-values.

        :param ax: matplotlib axis
        :param chrom: chromosome name
        :param region_start: start position of the region
        :param region_end: end position of the region
        :return: None
        """
        gwas_reader = ReadGwas(open(self.properties['file'], 'r'), has_header=self.properties['file_has_header'])

        # Fill in the position and pvalues lists with data from the GWAS file
        position = []
        pvalues = []
        for record in gwas_reader:
            if record.chromosome == chrom and region_start <= record.position <= region_end:
                position.append(record.position)
                pvalues.append(-np.log10(record.pvalue) if record.pvalue > 0 else 0)  # Notice the -log10 transformation

        # Plot the scatterplot
        ax.scatter(position, pvalues, s=self.properties['marker_size'], color=self.properties['color'], marker='o',
                   edgecolors='black', linewidths=.66)
