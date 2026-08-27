# Authors: Zepeng Mu (zmu@broadinstitute.org) and Yang I. Li (yangili1@uchicago.edu)
# Edits: Lucille Lopez-Delisle (lucille.delisle@unige.ch)

from .GenomeTrack import GenomeTrack, HUGE_NUMBER
import numpy as np
from ..utilities import plot_coverage, InputError, transform, change_chrom_names, opener, to_string, temp_file_from_intersect
import pyBigWig
from intervaltree import IntervalTree, Interval
import matplotlib
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from tqdm import tqdm
from . BigWigTrack import BigWigTrack

Path = mpath.Path

DEFAULT_LINKS_COLOR = 'blue'
DEFAULT_BIGWIG_COLOR = '#33a02c'


class SashimiBigWigTrack(BigWigTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = 'sashimiBigWig'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
bw_color = #666666
# To use a different color for negative values
#negative_color = red
# To use transparency, you can use alpha
# default is 1
#alpha = 0.5
# the default for min_value and max_value is 'auto' which means that the scale will go
# roughly from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# The number of bins takes the region to be plotted and divides it
# into the number of bins specified
# Then, at each bin the bigwig mean value is computed and plotted.
# A lower number of bins produces a coarser tracks
number_of_bins = 700
# to convert missing data (NaNs) into zeros. Otherwise, missing data is not plotted.
nans_to_zeros = true
# The possible summary methods are given by pyBigWig:
# mean/average/stdev/dev/max/min/cov/coverage/sum
# default is mean
summary_method = mean
# for type, the options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
#type = line:0.5
#type = points:0.5
# set show_data_range to false to hide the text on the left showing the data range
show_data_range = true
# To log transform your data you can also use transform and log_pseudocount:
# For the transform values:
# 'log1p': transformed_values = log(1 + initial_values)
# 'log': transformed_values = log(log_pseudocount + initial_values)
# 'log2': transformed_values = log2(log_pseudocount + initial_values)
# 'log10': transformed_values = log10(log_pseudocount + initial_values)
# '-log': transformed_values = - log(log_pseudocount + initial_values)
# For example:
#tranform = log
#log_pseudocount = 2
# When a transformation is applied, by default the y axis
# gives the transformed values, if you prefer to see
# the original values:
#y_axis_values = original
# If you want to have a grid on the y-axis
#grid = true
## Links customization
# The link file should be a BED file where the score
# is in the 5th column
#link_file =
# If the bed file contains a column for color (column 9), then this color can be used by
# setting:
#link_color = bed_rgb
# if link_color is a valid colormap name (like RbBlGn), then the score (column 5) is mapped
# to the colormap.
# In this case, the the link_min_value and link_max_value for the score can be provided, otherwise
# the maximum score and minimum score found are used.
#link_color = RdYlBu
#link_min_value=0
#link_max_value=100
# If the link_color is simply a color name, then this link_color is used and the score is not considered for the color.
link_color = darkblue
# To use transparency, you can use link_alpha
# default is 1
#link_alpha = 0.5
# options for link_line_style are 'solid', 'dashed', 'dotted', and 'dashdot'
link_line_style = solid
# The link in Sashimi plot is a Bezier curve.
# The height of the curve is calculated from the length of the intron.
# When the y-axis in bigwig track is different, the height of curve needs to be scaled.
link_scale_height = 1
# The line width for links is proportion to the numbers in the fifth column of the BED file (PSI).
# But the absolute width is calculated from the supplied numbers, which can look too thin or too wide sometimes.
# Use link_scale_line_width to scale the absolute line widths.
# You may need to try several values to get a satisfying result.
# Use this to scale Sashimi line width if the links are too thin or too wide.
#link_scale_line_width = 3
# Set link_line_width if you do not want width of links to scale with score (PSI).
# This overwrites link_scale_line_width.
#link_line_width = 2
# Set this to true to label scores (PSI) on links
link_labels = true
# optional: font size can be given to override the default size
#link_fontsize = 10
file_type = {TRACK_TYPE}
    """

    DEFAULTS_PROPERTIES = {
        # Bigwig related
        'max_value': None,
        'min_value': None,
        'show_data_range': True,
        'orientation': None,
        'bw_color': DEFAULT_BIGWIG_COLOR,
        'negative_color': None,
        'alpha': 1,
        'nans_to_zeros': False,
        'summary_method': 'mean',
        'number_of_bins': 700,
        'type': 'fill',
        'transform': 'no',
        'log_pseudocount': 0,
        'y_axis_values': 'transformed',
        'grid': False,
        # Links related
        'link_color': DEFAULT_LINKS_COLOR,
        'link_alpha': 1,
        'link_line_width': None,
        'link_line_style': 'solid',
        'link_max_value': None,
        'link_min_value': None,
        'link_scale_height': 1,
        'link_scale_line_width': 2,
        'link_labels': True,
        'link_fontsize': None,
        # General
        'region': None,  # Cannot be set manually but is set by tracksClass
    }
    NECESSARY_PROPERTIES = ['file', 'link_file']
    SYNONYMOUS_PROPERTIES = {
        'max_value': {
            'auto': None
        },
        'min_value': {
            'auto': None
        },
        'link_max_value': {
            'auto': None
        },
        'link_min_value': {
            'auto': None
        },
    }
    POSSIBLE_PROPERTIES = {
        'orientation': [None, 'inverted'],
        'summary_method': [
            'mean', 'average', 'max', 'min', 'stdev', 'dev', 'coverage', 'cov',
            'sum'
        ],
        'transform': ['no', 'log', 'log1p', '-log', 'log2', 'log10'],
        'y_axis_values': ['original', 'transformed'],
        'link_line_style': ['solid', 'dashed',
                            'dotted', 'dashdot'],
    }
    BOOLEAN_PROPERTIES = [
        'nans_to_zeros', 'show_data_range', 'grid', 'link_labels'
    ]
    STRING_PROPERTIES = [
        'file', 'file_type', 'overlay_previous', 'orientation',
        'summary_method', 'title', 'bw_color', 'negative_color', 'transform',
        'y_axis_values', 'type', 'link_file',
        'link_line_style', 'link_color'
    ]
    FLOAT_PROPERTIES = {
        'max_value': [-np.inf, np.inf],
        'min_value': [-np.inf, np.inf],
        'log_pseudocount': [-np.inf, np.inf],
        'alpha': [0, 1],
        'height': [0, np.inf],
        'link_fontsize': [0, np.inf],
        'link_alpha': [0, 1],
        'link_max_value': [-np.inf, np.inf],
        'link_min_value': [-np.inf, np.inf],
        'link_line_width': [0, np.inf],
        'link_scale_height': [0, np.inf],
        'link_scale_line_width': [0, np.inf]
    }
    INTEGER_PROPERTIES = {'number_of_bins': [1, np.inf]}

    # The bw_color can only be a color
    # negative_color can only be a color or None

    def __init__(self, *args, **kwargs):
        super(BigWigTrack, self).__init__(*args, **kwargs)
        self.bw = pyBigWig.open(self.properties['file'])
        self.show_number = self.properties['link_labels']

    def set_properties_defaults(self):
        super(BigWigTrack, self).set_properties_defaults()
        super(BigWigTrack, self).process_type_for_coverage_track()
        self.process_color('bw_color')
        if self.properties['negative_color'] is None:
            self.properties['negative_color'] = self.properties['bw_color']
        else:
            self.process_color('negative_color')
        # FOR LINK
        is_colormap = self.process_color('link_color', colormap_possible=True, bed_rgb_possible=True, default_value_is_colormap=False)
        self.interval_tree, min_score, max_score = self.process_bed(DEFAULT_LINKS_COLOR, file_key='link_file', color_key='link_color',
                                                                    plot_regions=self.properties['region'])
        # Initiate the colormap if needed
        self.colormap = None
        self.parametersUsingColormap = []
        if is_colormap:
            if self.properties['link_min_value'] is not None:
                min_score = self.properties['link_min_value']
            if self.properties['link_max_value'] is not None:
                max_score = self.properties['link_max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.colormap)
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
            self.parametersUsingColormap.append('link_color')

    def plot(self, ax, chrom_region, start_region, end_region):
        x_values, transformed_scores = self.get_transformed_values(chrom_region, start_region, end_region)
        if x_values is None:
            self.log.warning("Scores could not be computed. This will generate an empty track\n")
            return

        count = 0

        plot_coverage(ax, x_values, transformed_scores, self.plot_type,
                      self.size, self.properties['bw_color'],
                      self.properties['negative_color'],
                      self.properties['alpha'], self.properties['grid'])

        # Adjust ylim except for the inversion
        real_orientation = self.properties['orientation']
        self.properties['orientation'] = None

        self.adjust_ylim(ax)

        plot_ymin, plot_ymax = ax.get_ylim()
        self.pos_height = plot_ymax
        self.neg_height = plot_ymin

        # PLOT LINK
        if chrom_region not in self.interval_tree.keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in self.interval_tree.keys():
                self.log.warning("*Warning*\nNeither "
                                 + chrom_region_before + " nor "
                                 + chrom_region + " exists as a "
                                 "chromosome name inside the link_file."
                                 "No link will be plotted!!\n")
                self.interval_tree[chrom_region] = IntervalTree()
        arcs_in_region = sorted(
            self.interval_tree[chrom_region][start_region:end_region])
        for idx, interval in enumerate(arcs_in_region):
            # skip intervals whose start and end are outside the plotted region
            if interval.begin < start_region and interval.end > end_region:
                continue
            score_start = float(
                self.bw.values(chrom_region, interval.begin,
                               interval.begin + 1)[0])
            score_end = float(
                self.bw.values(chrom_region, interval.end,
                               interval.end + 1)[0])
            # Transform the scores
            score_start, score_end = \
                transform(np.array([score_start, score_end]),
                          self.properties['transform'],
                          self.properties['log_pseudocount'],
                          self.properties['file'])

            if self.properties['link_line_width'] is not None:
                self.line_width = float(self.properties['link_line_width'])
            else:
                self.line_width = self.properties['link_scale_line_width'] * np.log(
                    interval.data.score + 1) * 1.5

            self.plot_bezier(ax, interval, idx, score_start, score_end,
                             plot_ymin, plot_ymax)
            count += 1

        self.log.debug(f"{count} links plotted")

        # Adjust the ylim to include the potential links plotted
        # and use orientation
        if real_orientation == 'inverted':
            self.properties['orientation'] = 'inverted'
            ax.set_ylim(self.pos_height, self.neg_height)
        else:
            ax.set_ylim(self.neg_height, self.pos_height)

        return ax

    def plot_bezier(self, ax, interval, idx, start_height, end_height, ymin, ymax):

        def cubic_bezier(pts, t):
            b_x = (1 - t)**3 * pts[0][0] + 3 * t * (1 - t)**2 * pts[1][
                0] + 3 * t**2 * (1 - t) * pts[2][0] + t**3 * pts[3][0]
            b_y = (1 - t)**3 * pts[0][1] + 3 * t * (1 - t)**2 * pts[1][
                1] + 3 * t**2 * (1 - t) * pts[2][1] + t**3 * pts[3][1]
            return ((b_x, b_y))

        # width = (interval.end - interval.begin)

        height = (ymax - ymin) * 0.25 * self.properties['link_scale_height']
        epsilon = (ymax - ymin) * 0.05
        rgb = self.get_rgb(interval.data, param='link_color', default=DEFAULT_LINKS_COLOR)

        # Plot below x-axis
        if idx % 2 != 0:
            pts = [(interval.begin, ymin), (interval.begin, ymin-height),
                   (interval.end, ymin-height), (interval.end, ymin)]
            midpt = cubic_bezier(pts, 0.5)
            minpt = min(
                [cubic_bezier(pts, x)[1] for x in np.arange(0, 1, 0.05)])
            if minpt < self.neg_height:
                self.neg_height = minpt - epsilon

            pp1 = mpatches.PathPatch(Path(
                pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                     fc="none",
                                     ec=rgb,
                                     lw=self.line_width,
                                     ls=self.properties['link_line_style'],
                                     alpha=self.properties['link_alpha'])
            ax.add_patch(pp1)
            if self.show_number:
                ax.text(midpt[0],
                        midpt[1],
                        round(interval.data.score, 3),
                        fontsize=self.properties['link_fontsize'],
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=dict(facecolor='white', edgecolor='none', pad=0))
        # Plot above
        else:
            pts = [(interval.begin, start_height),
                   (interval.begin, height + start_height),
                   (interval.end, height + end_height),
                   (interval.end, end_height)]

            midpt = cubic_bezier(pts, 0.5)
            maxpt = max(
                [cubic_bezier(pts, x)[1] for x in np.arange(0, 1, 0.05)])
            if maxpt > self.pos_height:
                self.pos_height = maxpt + epsilon

            pp1 = mpatches.PathPatch(Path(
                pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                     fc="none",
                                     ec=rgb,
                                     lw=self.line_width,
                                     ls=self.properties['link_line_style'],
                                     alpha=self.properties['link_alpha'])
            ax.add_patch(pp1)
            if self.show_number:
                ax.text(midpt[0],
                        midpt[1],
                        round(interval.data.score, 3),
                        fontsize=self.properties['link_fontsize'],
                        horizontalalignment='center',
                        verticalalignment='center',
                        bbox=dict(facecolor='white', edgecolor='none', pad=0))

    # This y axis does not show the negative part, which is only Sashimi links
    def plot_y_axis(self,
                    ax,
                    plot_axis):
        """
        Plot the scale of the y axis with respect to the plot_axis
        Args:
            ax: axis to use to plot the scale
            plot_axis: the reference axis to get the max and min.

        Returns:

        """
        if not self.properties.get('show_data_range', True):
            return

        def value_to_str(value):
            # given a numeric value, returns a
            # string that removes unneeded decimal places
            if value % 1 == 0:
                str_value = str(int(value))
            else:
                str_value = f"{value:.1f}"
            return str_value

        def untransform(value, transform, log_pseudocount):
            # given a numeric value, transform and log_pseudocount
            # return the value before the transformation
            if transform == 'log':
                return np.exp(value) - log_pseudocount
            elif transform == 'log2':
                return np.exp2(value) - log_pseudocount
            elif transform == 'log10':
                return np.power(10, value) - log_pseudocount
            elif transform == 'log1p':
                return np.expm1(value)
            elif transform == '-log':
                return np.exp(-value) - log_pseudocount

        ymin, ymax = plot_axis.get_ylim()

        if not self.properties['grid']:
            if self.properties['min_value'] is None:
                min_value = 0
            else:
                min_value = self.properties['min_value']
            min_value_transformed = transform(np.array([min_value]),
                                            self.properties['transform'],
                                            self.properties['log_pseudocount'],
                                            self.properties['file'])[0]
            if self.properties['orientation'] is None:
                ymin = min_value_transformed
            else:
                ymax = min_value_transformed
        # If the ticks are closer than epsilon from the top or bottom
        # The vertical alignment of label is adjusted
        epsilon = (ymax - ymin) / 100
        # When the ymax and ymin are plotted (when there is no grid)
        # The tick is shifted inside of epsilon_pretty
        # To avoid to have only half of the width of the line plotted
        epsilon_pretty = epsilon
        y_axis = self.properties['y_axis_values']
        log_pseudocount = self.properties['log_pseudocount']
        transform_prop = self.properties['transform']

        if self.properties['grid']:
            # plot something that looks like this:
            # tick3 ┐
            #       │
            # tick2-|
            #       │
            # tick1 ┘
            if ymin < ymax:
                ticks_values = [
                    t for t in plot_axis.get_yticks()
                    if t <= ymax and t >= ymin
                ]
            else:
                ticks_values = [
                    t for t in plot_axis.get_yticks()
                    if t >= ymax and t <= ymin
                ]
                ticks_values.sort(reverse=True)
            labels_pos = ticks_values
            if transform_prop == 'no' or y_axis == 'transformed':
                ticks_labels = [value_to_str(t) for t in ticks_values]
            else:
                # There is a transformation and we want to display original values
                ticks_labels = [
                    value_to_str(untransform(t, transform_prop, log_pseudocount))
                    for t in ticks_values
                ]
        elif transform_prop == 'no' or y_axis == 'transformed':
            # This is a linear scale
            # plot something that looks like this:
            # ymax ┐
            #      │
            #      │
            # ymin ┘
            # adjust the positions such that the lines are plotted complete
            # and not only half of the width of the line.
            ticks_values = [ymin + epsilon_pretty, ymax - epsilon_pretty]
            labels_pos = [ymin, ymax]
            ticks_labels = [value_to_str(v) for v in [ymin, ymax]]
            if y_axis == 'transformed' and transform_prop != 'no':
                if transform_prop == 'log1p':
                    ymid_str = "log(1 + x)"
                else:
                    if log_pseudocount == 0:
                        ymid_str = f"{transform_prop}(x)"
                    else:
                        ymid_str = f"{transform_prop}({log_pseudocount} + x)"

                ax.text(0, (ymax + ymin) / 2,
                        ymid_str,
                        verticalalignment='center',
                        horizontalalignment='right',
                        wrap=True)
        else:
            # There is a transformation and we want to display original values
            if ymin * ymax < 0:
                ymid = 0
            else:
                ymid = (ymin + ymax) / 2
            # plot something that looks like this:
            # ymax ┐
            #      │
            # ymid-|
            #      │
            # ymin ┘
            ticks_values = [ymin + epsilon_pretty, ymid, ymax - epsilon_pretty]
            labels_pos = [ymin, ymid, ymax]
            ticks_labels = [
                value_to_str(untransform(v, transform_prop, log_pseudocount))
                for v in [ymin, ymid, ymax]
            ]

        # The lower label should be verticalalignment='bottom'
        # if it corresponds to ymin
        i = 0
        if (ymin < ymax and ticks_values[i] <= ymin + epsilon) \
           or (ymin > ymax and ticks_values[i] >= ymin + epsilon):
            v_al = 'bottom'
            adjusted_value = labels_pos[i] - epsilon
        else:
            v_al = 'center'
            adjusted_value = labels_pos[i]
        ax.text(-0.2,
                adjusted_value,
                ticks_labels[i],
                verticalalignment=v_al,
                horizontalalignment='right')
        x_pos = [0, 0.5]
        y_pos = [ticks_values[i]] * 2
        for i in range(1, len(ticks_values) - 1):
            ax.text(-0.2,
                    labels_pos[i],
                    ticks_labels[i],
                    verticalalignment='center',
                    horizontalalignment='right')
            x_pos += [0.5, 0, 0.5]
            y_pos += [ticks_values[i]] * 3

        # The upper label should be verticalalignment='top'
        # if it corresponds to ymax
        i = len(ticks_values) - 1
        if (ymin < ymax and ticks_values[i] >= ymax - epsilon) \
           or (ymin > ymax and ticks_values[i] <= ymax - epsilon):
            v_al = 'top'
        else:
            v_al = 'center'
        ax.text(-0.2,
                labels_pos[i],
                ticks_labels[i],
                verticalalignment=v_al,
                horizontalalignment='right')
        x_pos += [0.5, 0]
        y_pos += [ticks_values[i]] * 2

        # Finally plot the line:
        ax.plot(x_pos, y_pos, color='black', linewidth=1)

        # Set the lims:
        ax.set_ylim(plot_axis.get_ylim())
        ax.set_xlim(0, 1)
        ax.patch.set_visible(False)

    def get_scores(self, bw_var, file, chrom_region, start_region, end_region):
        bw = eval(bw_var)
        scores_per_bin = None
        if chrom_region not in bw.chroms().keys():
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in bw.chroms().keys():
                self.log.warning("*Warning*\nNeither " + chrom_region_before +
                                 " nor " + chrom_region + " exists as a "
                                 "chromosome name inside the bigwig file. "
                                 "No score will be computed for"
                                 f" {file}.\n")
                scores_per_bin = np.array([np.nan] *
                                          self.properties['number_of_bins'])

        if scores_per_bin is None and start_region > bw.chroms()[chrom_region]:
            self.log.warning(
                "*Warning*\nThe region to plot starts beyond the"
                " chromosome size. No score will be computed for"
                f" {file}.\n"
                f"{chrom_region} size: {bw.chroms()[chrom_region]}"
                f". Region to plot {start_region}-{end_region}\n")
            scores_per_bin = np.array([np.nan] *
                                      self.properties['number_of_bins'])

        if scores_per_bin is None and end_region > bw.chroms()[chrom_region]:
            self.log.warning(
                "*Warning*\nThe region to plot extends beyond the"
                " chromosome size. Please check.\n"
                f"{chrom_region} size: {bw.chroms()[chrom_region]}"
                f". Region to plot {start_region}-{end_region}\n")
            temp_end_region = bw.chroms()[chrom_region]
            temp_nbins = int(self.properties['number_of_bins'] *
                             (temp_end_region - start_region) /
                             (end_region - start_region))
        else:
            temp_end_region = end_region
            temp_nbins = self.properties['number_of_bins']
        # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        # of the memory. This only occurs when calling trackPlot from different
        # processors. Reloading the file solves the problem.
        if scores_per_bin is None:
            num_tries = 0
            while num_tries < 5:
                num_tries += 1
                try:
                    scores_per_bin = np.array(
                        bw.stats(
                            chrom_region,
                            start_region,
                            temp_end_region,
                            nBins=temp_nbins,
                            type=self.properties['summary_method'])).astype(
                                float)
                except Exception as e:
                    bw = pyBigWig.open(self.properties['file'])

                    self.log.warning("error found while reading bigwig scores "
                                     f"({e}).\nTrying again."
                                     f" Iter num: {num_tries}.\n")
                    pass
                else:
                    if num_tries > 1:
                        self.log.warning(
                            f"After {num_tries} the scores could be computed.\n"
                        )
                    break
        return temp_end_region, temp_nbins, scores_per_bin

    def __del__(self):
        try:
            self.bw.close()
        except AttributeError:
            pass
