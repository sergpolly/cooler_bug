from cooler import Cooler

# input_uri = "unzoomifiable_5kb.cool"
input_uri = "missing_bin_cooler_40kb.cool"

clr = Cooler(input_uri)

print(clr.pixels(join=True)[17872226])


table = clr.pixels(join=True)[17800000:]

print(table[table['chrom2'].isnull()])
#          chrom1   start1     end1 chrom2  start2  end2  count
# 17872226  chr88  2640000  2680000    NaN     NaN   NaN      1


########################
# TODO
#########################
# Now i'm on the same page with Nezar
# bad pixel is found ...
# FIND OUT  - HOW that happened ?!?!?!?!?! ...
# 
# processes that might be involved:
# 
#     pairix ${pairs_lib}
# 
#     cooler cload pairix \
#        --nproc ${task.cpus} \
#        --assembly ${params.input.genome.assembly} \
#        ${chrom_sizes}:${res} ${pairs_lib} ${library}.${res}.cool
