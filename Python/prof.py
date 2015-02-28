import pstats
p=pstats.Stats('out.tmp')
p.strip_dirs().sort_stats('cumulative').print_stats(20)