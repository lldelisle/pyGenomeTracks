from . GenomeTrack import GenomeTrack
from pygenometracks.utilities import InputError
from intervaltree import IntervalTree, Interval
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon
from .. utilities import opener, to_string, change_chrom_names, temp_file_from_intersect
from tqdm import tqdm

DEFAULT_LINKS_COLOR = 'blue'
HUGE_NUMBER = 1e15  # Which should be above any chromosome size


class LinksTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.arcs', '.arc', '.link', '.links', '.bedpe']
    TRACK_TYPE = 'links'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + f"""
# the file format for links is (tab separated)
#   chr1 start1 end1 chr2 start2 end2 (score ...)
# The score field is optional
# The fields after the score field will be ignored
# for example:
#   chr1 100 200 chr1 250 300 0.5
# depending on the value of links_type either 'arcs' or 'triangles' or 'loops' can be plotted.
# If arcs, a line will be drawn from the center of the first region (chr1: 150),
# to the center of the other region (chr1: 275).
# if triangles, the vertix of the triangle will be drawn at the center between the two points
# (also the extremity of each position is used)
# if loops, a rectangle highlighting the intersection between the 2 regions will be shown
# the triangles, and loops options are convenient to overlay over a
# Hi-C matrix to highlight the matrix pixel of the highlighted link
# For these tracks do not hesitate to put large line_width like 5 or 10.
links_type = arcs
# For triangles and arcs, by default the extremities coordinates are used
# To use the middle of start1 and end1 and the middle of start2 and end2
#use_middle = true
# color of the lines
color = red
# if color is a valid colormap name (like RdYlGn),
# then the score is mapped to the colormap.
#color = RdYlGn
# To set the minimum and maximum value of the colormap:
#min_value = 0
#max_value = 1.2
# To use transparency, you can use alpha
# default is 0.8
# alpha = 0.5
# if line_width is not given, the score is used to set the line width
# using the following formula (0.5 * square root(score)
#line_width = 0.5
# options for line_style are 'solid', 'dashed', 'dotted', and 'dashdot'
line_style = solid
# If you want to compact the arcs (when you have both long and short arcs)
# You can choose a compact level of
# 1 (the height is proportional to the square root of the distance)
# 2 (the height is the same for all distances)
# (default is 0 proportional to distance)
#compact_arcs_level = 2
# To be able to see small arcs when big arcs exists, you can set
# the upper y limit.
# The unit is bp. This corresponds to the longest arc you will see.
# This option is incompatible with compact_arcs_level = 2
#ylim = 100000
file_type = {TRACK_TYPE}
    """
    DEFAULTS_PROPERTIES = {'links_type': 'arcs',
                           'line_width': None,
                           'line_style': 'solid',
                           'orientation': None,
                           'color': DEFAULT_LINKS_COLOR,
                           'alpha': 0.8,
                           'max_value': None,
                           'min_value': None,
                           'region': None,  # Cannot be set manually but is set by tracksClass
                           'ylim': None,
                           'compact_arcs_level': '0',
                           'use_middle': False}
    NECESSARY_PROPERTIES = ['file']
    SYNONYMOUS_PROPERTIES = {'max_value': {'auto': None},
                             'min_value': {'auto': None},
                             'ylim': {'auto': None}}
    POSSIBLE_PROPERTIES = {'orientation': [None, 'inverted'],
                           'links_type': ['arcs', 'triangles', 'loops'],
                           'line_style': ['solid', 'dashed',
                                          'dotted', 'dashdot'],
                           'compact_arcs_level': ['0', '1', '2']}
    BOOLEAN_PROPERTIES = ['use_middle']
    STRING_PROPERTIES = ['file', 'file_type', 'overlay_previous',
                         'orientation', 'links_type', 'line_style',
                         'title', 'color', 'compact_arcs_level']
    FLOAT_PROPERTIES = {'max_value': [- np.inf, np.inf],
                        'min_value': [- np.inf, np.inf],
                        'ylim': [0, np.inf],
                        'alpha': [0, 1],
                        'line_width': [0, np.inf],
                        'height': [0, np.inf]}
    INTEGER_PROPERTIES = {}
    # The color can be a color or a colormap (if there is a score)

    def set_properties_defaults(self):
        super(LinksTrack, self).set_properties_defaults()
        self.max_height = None
        self.interval_tree, min_score, max_score, has_score = self.process_link_file(self.properties['region'])
        if self.properties['line_width'] is None and not has_score:
            self.log.warning("*WARNING* for section "
                             f"{self.properties['section_name']}"
                             " no line_width has been set but some "
                             "lines do not have scores."
                             "line_width has been set to "
                             "0.5.\n")
            self.properties['line_width'] = 0.5

        self.colormap = None
        # check if the color given is a color map
        is_colormap = self.process_color('color', colormap_possible=True,
                                         default_value_is_colormap=False)
        if is_colormap:
            if not has_score:
                self.log.warning("*WARNING* for section "
                                 f"{self.properties['section_name']}"
                                 " a colormap was chosen but some "
                                 "lines do not have scores."
                                 "Color has been set to "
                                 f"{DEFAULT_LINKS_COLOR}.\n")
                self.properties['color'] = DEFAULT_LINKS_COLOR
            else:
                self.colormap = self.properties['color']

        if self.colormap is not None:
            if self.properties['min_value'] is not None:
                min_score = self.properties['min_value']
            if self.properties['max_value'] is not None:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

        if self.properties['compact_arcs_level'] == '2' and \
           self.properties['ylim'] is not None:
            self.log.warning("*WARNING* for section "
                             f"{self.properties['section_name']}"
                             " a ylim was set but "
                             "compact_arcs_level was set to 2."
                             "ylim will be ignore.\n")
            self.properties['ylim'] = None

    def plot(self, ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        :param ax: matplotlib axis
        :param label_ax: matplotlib axis for labels
        """
        self.max_height = 0
        count = 0

        if chrom_region not in list(self.interval_tree):
            chrom_region_before = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            if chrom_region not in list(self.interval_tree):
                self.log.warning("*Warning*\nNeither " + chrom_region_before
                                 + " nor " + chrom_region + " exists as a "
                                 "chromosome name inside the link file. "
                                 "This will generate an empty track!!\n")
                return

        chrom_region = self.check_chrom_str_bytes(self.interval_tree, chrom_region)

        arcs_in_region = sorted(self.interval_tree[chrom_region][region_start:region_end])

        for idx, interval in enumerate(arcs_in_region):
            # skip intervals whose start and end are outside the plotted region
            if interval.begin < region_start and interval.end > region_end:
                continue

            if self.properties['line_width'] is not None:
                self.line_width = float(self.properties['line_width'])
            else:
                self.line_width = 0.5 * np.sqrt(interval.data[4])

            if self.properties['links_type'] == 'triangles':
                self.plot_triangles(ax, interval)
            elif self.properties['links_type'] == 'loops':
                self.plot_loops(ax, interval.data)
            else:
                self.plot_arcs(ax, interval)

            count += 1

        # the arc height is equal to the radius, the track height is the largest
        # radius plotted plus an small increase to avoid cropping of the arcs
        self.max_height *= 1.1
        if self.properties['ylim'] is None:
            ymax = self.max_height
        else:
            if self.properties['compact_arcs_level'] == '1':
                ymax = np.sqrt(self.properties['ylim'])
            else:
                ymax = self.properties['ylim']
        self.log.debug(f"{count} were links plotted")
        if self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, -1)
        else:
            ax.set_ylim(-1, ymax)

    def plot_y_axis(self, ax, plot_ax):
        if self.colormap is not None and self.properties['overlay_previous'] == 'no':
            self.colormap.set_array([])

            cobar = plt.colorbar(self.colormap, ax=ax, fraction=1,
                                 orientation='vertical')

            cobar.solids.set_edgecolor("face")
            cobar.ax.tick_params(labelsize='smaller')
            cobar.ax.yaxis.set_ticks_position('left')
            # adjust the labels of the colorbar
            ticks = cobar.ax.get_yticks()
            labels = cobar.ax.set_yticklabels(ticks.astype('float32'))
            (vmin, vmax) = cobar.mappable.get_clim()
            for idx in np.where(ticks == vmin)[0]:
                # if the label is at the start of the colobar
                # move it above avoid being cut or overlapping with other track
                labels[idx].set_verticalalignment('bottom')
            for idx in np.where(ticks == vmax)[0]:
                # if the label is at the end of the colobar
                # move it a bit inside to avoid overlapping
                # with other labels
                labels[idx].set_verticalalignment('top')

    def plot_arcs(self, ax, interval):

        width = (interval.end - interval.begin)
        if self.properties['compact_arcs_level'] == '1':
            half_height = np.sqrt(width)
        elif self.properties['compact_arcs_level'] == '2':
            half_height = 1000
        else:
            half_height = width
        center = interval.begin + width / 2
        if half_height > self.max_height:
            self.max_height = half_height
        # I think this was an error...
        # ax.plot([center], [diameter])
        if self.colormap:
            # translate score field
            # into a color
            rgb = self.colormap.to_rgba(interval.data[4])
        else:
            rgb = self.properties['color']
        ax.add_patch(Arc((center, 0), width,
                         2 * half_height, 0, 0, 180, color=rgb,
                         linewidth=self.line_width,
                         ls=self.properties['line_style']))

    def plot_triangles(self, ax, interval):
        x1 = interval.begin
        x2 = x1 + float(interval.end - interval.begin) / 2
        x3 = interval.end
        y1 = 0
        if self.properties['compact_arcs_level'] == '1':
            y2 = np.sqrt(interval.end - interval.begin)
        elif self.properties['compact_arcs_level'] == '2':
            y2 = 1000
        else:
            y2 = (interval.end - interval.begin)

        if self.colormap:
            # translate score field
            # into a color
            rgb = self.colormap.to_rgba(interval.data[4])
        else:
            rgb = self.properties['color']

        triangle = Polygon(np.array([[x1, y1], [x2, y2], [x3, y1]]),
                           closed=False,
                           facecolor='none', edgecolor=rgb,
                           linewidth=self.line_width,
                           ls=self.properties['line_style'])
        ax.add_artist(triangle)
        if y2 > self.max_height:
            self.max_height = y2

    def plot_loops(self, ax, loop):
        """
              " <- 2
        3->  "  " <- 1
               " <- 0
            """
        width1 = loop[1] - loop[0]
        width2 = loop[3] - loop[2]
        x0 = (loop[1] + loop[2]) / 2
        y0 = loop[2] - loop[1]

        x1 = x0 + width2 / 2
        y1 = y0 + width2

        x2 = (loop[0] + loop[3]) / 2
        y2 = loop[3] - loop[0]

        x3 = x0 - width1 / 2
        y3 = y0 + width1

        if self.colormap:
            # translate score field
            # into a color
            rgb = self.colormap.to_rgba(loop[4])
        else:
            rgb = self.properties['color']

        rectangle = Polygon(np.array([[x0, y0], [x1, y1], [x2, y2], [x3, y3]]),
                            facecolor='none', edgecolor=rgb,
                            linewidth=self.line_width,
                            ls=self.properties['line_style'])
        ax.add_artist(rectangle)
        if y2 > self.max_height:
            self.max_height = y2

    def process_link_file(self, plot_regions):
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is a score.

        if plot_regions is None:
            file_to_open = self.properties['file']
        else:
            # To be sure we do not miss links we will intersect with bed with
            # only chromosomes used in plot_regions
            plot_regions_adapted = [(chrom, 0, HUGE_NUMBER) for chrom, __, __ in plot_regions]
            file_to_open = temp_file_from_intersect(self.properties['file'],
                                                    plot_regions_adapted)

        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        has_score = True
        max_score = float('-inf')
        min_score = float('inf')
        file_h = opener(file_to_open)
        for line in tqdm(file_h.readlines()):
            line_number += 1
            line = to_string(line)
            if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                continue
            try:
                chrom1, start1, end1, chrom2, start2, end2 = line.strip().split('\t')[:6]
            except Exception as detail:
                raise InputError('File not valid. The format is chrom1'
                                 ' start1, end1, '
                                 f'chrom2, start2, end2\nError: {detail}\n'
                                 f' in line\n {line}')
            if chrom1 != chrom2:
                self.log.warning(f"Only links in same chromosome are used. Skipping line\n{line}\n")
                continue

            try:
                score = line.strip().split('\t')[6]
            except IndexError:
                has_score = False
                score = np.nan

            try:
                start1 = int(start1)
                end1 = int(end1)
                start2 = int(start2)
                end2 = int(end2)
            except ValueError as detail:
                raise InputError(f"Error reading line: {line_number}. One of the fields is not "
                                 f"an integer.\nError message: {detail}")

            assert start1 <= end1, f"Error in line #{line_number}, end1 larger than start1 in {line}"
            assert start2 <= end2, f"Error in line #{line_number}, end2 larger than start2 in {line}"

            if has_score:
                try:
                    score = float(score)
                except ValueError as detail:
                    self.log.warning(f"Warning: reading line: {line}. The score is not valid {score} will not be used. "
                                     f"\nError message: {detail}\n")
                    score = np.nan
                    has_score = False
                else:
                    if score < min_score:
                        min_score = score
                    if score > max_score:
                        max_score = score

            if chrom1 not in interval_tree:
                interval_tree[chrom1] = IntervalTree()

            if start2 < start1:
                start1, start2 = start2, start1
                end1, end2 = end2, end1

            if self.properties['use_middle']:
                mid1 = (start1 + end1) / 2
                mid2 = (start2 + end2) / 2
                interval_tree[chrom1].add(Interval(mid1, mid2, [start1, end1, start2, end2, score]))
            else:
                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, [start1, end1, start2, end2, score]))
            valid_intervals += 1

        if valid_intervals == 0:
            self.log.warning(f"No valid intervals were found in file {self.properties['file']}.\n")

        file_h.close()
        return(interval_tree, min_score, max_score, has_score)
