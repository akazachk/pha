from lgreader import LGReader; 
reader = LGReader("/Users/akazachk/repos/pha2/results/saved_results/pha.csv")
reader.set_param("CUT_LIMIT",1000);
reader.set_param("CUT_PRESOLVE",0);
reader.set_inst("stein15_nosym", compute_gap_closed=False)
reader.hplane_analysis()
