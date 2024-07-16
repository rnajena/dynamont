awk '{print $0 "\t" (NR==1 ? "level_stdv" :"0.5")}' rna_r9.4_180mv_70bps_extended.model > rna_r9.4_180mv_70bps_extended_stdev0_5.model
awk '{print $0 "\t" (NR==1 ? "level_stdv" :"1.0")}' rna_r9.4_180mv_70bps_extended.model > rna_r9.4_180mv_70bps_extended_stdev1.model
