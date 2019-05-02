import numpy as np                    
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle
def draw_gene_simple(ax, bed, ypos, patch_height = None, arrow_length= None,
                     rgb=None, edgecolor=None, linewidth=None, 
                     **kwargs
                    ):
    """
    draws an interval with direction (if given)
    """
    if arrow_length is None:
        arrow_length = np.diff(ax.get_xlim())[0]/20.
    if patch_height is None:
        patch_height = np.diff(ax.get_ylim())[0]/5.
        
    from matplotlib.patches import Polygon

    if bed.strand not in ['+', '-']:
        ax.add_patch(Rectangle((bed.start, ypos), bed.end - bed.start, patch_height,
                               edgecolor=edgecolor, facecolor=rgb, linewidth=linewidth,
                               **kwargs))
    else:

        vertices = _draw_arrow(bed.start, bed.end, bed.strand, ypos, patch_height, arrow_length )
        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=rgb,
                             linewidth=linewidth,
                             **kwargs))
#     _draw_arrow

def _draw_arrow(start, end, strand, ypos, height, arrow_length):
    """
    Draws a filled arrow
    :param start:
    :param end:
    :param strand:
    :param ypos:
    :return: None
    """
    half_height = float(height) / 2
    if strand == '+':
        x0 = start
        x1 = end  # - self.small_relative
        y0 = ypos
        y1 = ypos + height
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        -----------------\
        ---------------- /

        """

        vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + arrow_length, y0 + half_height), (x1, y0)]

    else:
        x0 = start  # + self.small_relative
        x1 = end
        y0 = ypos
        y1 = ypos + height
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        /-----------------
        \-----------------
        """
        vertices = [(x0, y0), (x0 - arrow_length, y0 + half_height), (x0, y1), (x1, y1), (x1, y0)]

    return vertices

