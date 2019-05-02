# %pdb 0
import sys,os
os.chdir(__file__.rsplit('/',1)[0])
# print os.listdir('.')
import pygenometracks.tracksClass
import pygenometracks as pygtk
import pygenometracks.tracks.BedTrack as btk
# btk = pygtk.tracks.BedTrack.__module__
import pygenometracks.plotters as plotters
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2,1,sharex='col')

# import 
# tk = pygtk.tracksClass.BigWigTrack(dict(file='test-data/0308__bam-to-bw/ALIGNMENT-sorted.bw'))
# tk = pygtk.tracksClass.DevBigWigTrack(dict(file='test-data/0308__bam-to-bw/ALIGNMENT-sorted.bw'))
# tk.plot(axs[0],*region)


region = 'Chr1',13898009, 13899205
chrom_region,start_region,end_region = region
ax = axs[1]

def getBedIntervals():
    self = tk_bed = tk = pygtk.tracksClass.BedTrack( 
        dict(file = 'test-plotter.bed6',section_name='wtf',title='title',interval_height=0.03 ) ,)
    self.small_relative = 0.004 * (end_region - start_region)
    self.properties['interval_height'] = 0.05
    tk.process_bed()
    genes_overlap = tk.interval_tree[region[0]][region[1]:region[2]]
    genes_overlap = sorted(genes_overlap)
    genes_overlap = [x.data for x in genes_overlap]
#     d = genes_overlap[0].data
    return genes_overlap

genes_overlap = getBedIntervals()

class temp(object):
    pass
d = temp()
vars(d).update(vars(genes_overlap[0]))
# d = pyutil.util_obj(**vars(d))
iv = d



ax.set_xlim(*region[1:])

iv.strand='.'
plotters.draw_gene_simple(ax,  iv, 0.75, )

iv.strand='+'
plotters.draw_gene_simple(ax,  iv, 0.25, )

fig.savefig(__file__.rsplit('/',1)[1] +'.png')