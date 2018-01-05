from cooler.io import ContactBinner
from cooler.util import binnify, get_binsize, parse_region
from cooler import Cooler
from cooler.tools import lock


import itertools
import six
import numpy as np
import pandas



def check_bins(bins, chromsizes):
    is_cat = pandas.api.types.is_categorical(bins['chrom'])
    bins = bins.copy()
    if not is_cat:
        bins['chrom'] = pandas.Categorical(
            bins.chrom, 
            categories=list(chromsizes.index), 
            ordered=True)
    else:
        assert (bins['chrom'].cat.categories == chromsizes.index).all()

    return bins



class GenomeSegmentation(object):
    def __init__(self, chromsizes, bins):
        bins = check_bins(bins, chromsizes)
        self._bins_grouped = bins.groupby('chrom', sort=False)
        nbins_per_chrom = self._bins_grouped.size().values

        self.chromsizes = chromsizes
        self.binsize = get_binsize(bins)
        self.contigs = list(chromsizes.keys())
        self.bins = bins
        self.idmap = pandas.Series(
            index=chromsizes.keys(), 
            data=range(len(chromsizes)))
        self.chrom_binoffset = np.r_[0, np.cumsum(nbins_per_chrom)]
        self.chrom_abspos = np.r_[0, np.cumsum(chromsizes.values)]
        self.start_abspos = (self.chrom_abspos[bins['chrom'].cat.codes] + 
                             bins['start'].values)
    
    def fetch(self, region):
        chrom, start, end = parse_region(region, self.chromsizes)
        result = self._bins_grouped.get_group(chrom)
        if start > 0 or end < self.chromsizes[chrom]:
            lo = result['end'].values.searchsorted(start, side='right')
            hi = lo + result['start'].values[lo:].searchsorted(end, side='left')
            result = result.iloc[lo:hi]
        return result



class CoolerAggregator(ContactBinner):
    """
    Aggregate contacts from an existing Cooler file.
    """
    def __init__(self, source_uri, bins, chunksize, batchsize, map=map):
        from cooler.api import Cooler
        self._map = map
        self.source_uri = source_uri
        self.chunksize = chunksize
        self.batchsize = batchsize

        clr = Cooler(source_uri)
        self._size = clr.info['nnz']
        self.old_binsize = clr.binsize
        self.old_chrom_offset = clr._load_dset('indexes/chrom_offset')
        self.old_bin1_offset = clr._load_dset('indexes/bin1_offset')
        self.gs = GenomeSegmentation(clr.chromsizes, bins)
        self.new_binsize = get_binsize(bins)
        assert self.new_binsize % self.old_binsize == 0
        self.factor = self.new_binsize // self.old_binsize
    
    def _aggregate(self, span):
        from cooler.api import Cooler
        lo, hi = span

        clr = Cooler(self.source_uri)
        # convert_enum=False returns chroms as raw ints
        table = clr.pixels(join=True, convert_enum=False)
        chunk = table[lo:hi]
        # logger.info('{} {}'.format(lo, hi))
        print('{} {}'.format(lo, hi))

        # use the "start" point as anchor for re-binning
        # XXX - alternatives: midpoint anchor, proportional re-binning
        binsize = self.gs.binsize
        chrom_binoffset = self.gs.chrom_binoffset
        chrom_abspos = self.gs.chrom_abspos
        start_abspos = self.gs.start_abspos

        chrom_id1 = chunk['chrom1'].values
        chrom_id2 = chunk['chrom2'].values
        start1 = chunk['start1'].values
        start2 = chunk['start2'].values
        if binsize is None:
            abs_start1 = chrom_abspos[chrom_id1] + start1
            abs_start2 = chrom_abspos[chrom_id2] + start2
            chunk['bin1_id'] = np.searchsorted(
                start_abspos, 
                abs_start1, 
                side='right') - 1
            chunk['bin2_id'] = np.searchsorted(
                start_abspos, 
                abs_start2, 
                side='right') - 1
        else:
            rel_bin1 = np.floor(start1/binsize).astype(int)
            rel_bin2 = np.floor(start2/binsize).astype(int)
            chunk['bin1_id'] = chrom_binoffset[chrom_id1] + rel_bin1
            chunk['bin2_id'] = chrom_binoffset[chrom_id2] + rel_bin2

        grouped = chunk.groupby(['bin1_id', 'bin2_id'], sort=False)
        return grouped['count'].sum().reset_index()

    def aggregate(self, span):
        try:
            chunk = self._aggregate(span)
        except MemoryError as e:
            raise RuntimeError(str(e))
        return chunk

    def __iter__(self):
        old_chrom_offset = self.old_chrom_offset
        old_bin1_offset = self.old_bin1_offset
        chunksize = self.chunksize
        batchsize = self.batchsize
        factor = self.factor
        
        # Partition pixels into chunks, respecting chrom1 boundaries
        spans = []
        for chrom, i in six.iteritems(self.gs.idmap):
            # it's important to extract some multiple of `factor` rows at a time
            c0 = old_chrom_offset[i]
            c1 = old_chrom_offset[i + 1]
            step = (chunksize // factor) * factor
            edges = np.arange(
                old_bin1_offset[c0], 
                old_bin1_offset[c1] + step, 
                step)
            edges[-1] = old_bin1_offset[c1]
            spans.append(zip(edges[:-1], edges[1:]))
        spans = list(itertools.chain.from_iterable(spans))
        
        # Process batches of k chunks at a time, then yield the results
        for i in range(0, len(spans), batchsize):
            try:
                lock.acquire()
                print("right before collapse ...{}  {}".format(i,spans[i:i+batchsize]))
                results = self._map(self.aggregate, spans[i:i+batchsize])
            finally:
                lock.release()
            for df in results:
                # yield {k: v.values for k, v in six.iteritems(df)}
                yield df





input_uri = ""

c = Cooler(input_uri)

new_bins = binnify(
                c.chromsizes,
                2*c.binsize)



iterator = CoolerAggregator(
                    input_uri,
                    new_bins,
                    1000000,
                    batchsize=1,
                    map=map)

# # last message before it fails ...
# # INFO:cooler:17868809 17872380
# for ii in iterator:
#     print(ii)




# from cooler.api import Cooler
lo, hi = 17869999, 17872300
# lo, hi = 17868809, 17872380

clr = Cooler(input_uri)
# convert_enum=False returns chroms as raw ints
table = clr.pixels(join=True, convert_enum=False)
chunk = table[lo:hi]
# logger.info('{} {}'.format(lo, hi))
print('{} {}'.format(lo, hi))

# use the "start" point as anchor for re-binning
# XXX - alternatives: midpoint anchor, proportional re-binning
binsize = iterator.gs.binsize
chrom_binoffset = iterator.gs.chrom_binoffset
chrom_abspos = iterator.gs.chrom_abspos
start_abspos = iterator.gs.start_abspos

chrom_id1 = chunk['chrom1'].values
chrom_id2 = chunk['chrom2'].values
start1 = chunk['start1'].values
start2 = chunk['start2'].values
if binsize is None:
    abs_start1 = chrom_abspos[chrom_id1] + start1
    abs_start2 = chrom_abspos[chrom_id2] + start2
    chunk['bin1_id'] = np.searchsorted(
        start_abspos, 
        abs_start1, 
        side='right') - 1
    chunk['bin2_id'] = np.searchsorted(
        start_abspos, 
        abs_start2, 
        side='right') - 1
else:
    rel_bin1 = np.floor(start1/binsize).astype(int)
    rel_bin2 = np.floor(start2/binsize).astype(int)
    chunk['bin1_id'] = chrom_binoffset[chrom_id1] + rel_bin1
    chunk['bin2_id'] = chrom_binoffset[chrom_id2] + rel_bin2

grouped = chunk.groupby(['bin1_id', 'bin2_id'], sort=False)
return grouped['count'].sum().reset_index()











