"""
This python script will generate two files:
- docs/content/all_default_properties_rst.txt
This file is a rst table with all the defaults values for each parameter
for each track class. This table is included in the readthedocs
- docs/content/all_possible_properties.txt
This file is a markdown list with possible values.
This can also be used in the readthedocs
"""
from pygenometracks.tracksClass import PlotTracks, DEFAULT_TRACK_HEIGHT
import numpy as np
import os.path

not_used_string = ''
used_string = 'X'
track_separator = ', '
defaut_for_all = {'overlay_previous': 'no'}


def main():
    all_tracks = PlotTracks.get_available_tracks()
    my_prefered_order_tracks_names = ['x_axis', 'epilogos', 'links',
                                      'domains', 'bed', 'gtf', 'narrow_peak',
                                      'bigwig', 'bedgraph', 'bedgraph_matrix',
                                      'hlines', 'hic_matrix', 'scalebar']
    my_prefered_order_tracks_names = [k for k in my_prefered_order_tracks_names
                                      if k in all_tracks]
    other_tracks = list(set(all_tracks.keys())
                        - set(my_prefered_order_tracks_names))
    # Get all possible and default parameters
    all_default_parameters = {}
    all_tracks_with_default = []
    for track_type in my_prefered_order_tracks_names + other_tracks:
        track_class = all_tracks[track_type]
        for p, value in defaut_for_all.items():
            all_default_parameters[p] = all_default_parameters.get(p, {})
            all_default_parameters[p][track_type] = value
        has_default = False
        for p, value in track_class.DEFAULTS_PROPERTIES.items():
            if p != 'region':
                all_default_parameters[p] = all_default_parameters.get(p, {})
                all_default_parameters[p][track_type] = value
                has_default = True
        if has_default:
            all_tracks_with_default += [track_type]

    # For the default they are summarized in a matrix
    mat = np.empty((len(all_default_parameters) + 1, len(all_tracks_with_default) + 1),
                   dtype='U25')
    mat[0, 0] = 'parameter'
    for j, track_type in enumerate(all_tracks_with_default, start=1):
        if j == mat.shape[1] - 1:
            mat[0, j] = track_type.replace('_', '\_') + ' \\\\'
        else:
            mat[0, j] = track_type.replace('_', '\_')

        for i, p in enumerate(all_default_parameters):
            if j == 1:
                mat[i + 1, 0] = p.replace('_', '\_')

            if track_type in all_default_parameters[p]:
                value = used_string
            else:
                value = not_used_string
            if j == mat.shape[1] - 1:
                value = value + ' \\\\'

            mat[i + 1, j] = value
    # The matrix is written in a file to be able to use it in latex
    # max_char = max([len(mat[i, j]) for i in range(mat.shape[0]) for j in range(mat.shape[1])])
    np.savetxt(os.path.join("param_per_track.txt"),
               mat, fmt=f'%s', delimiter=" & ",
               comments='')


if __name__ == "__main__":
    main()
