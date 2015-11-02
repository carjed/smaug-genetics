main_rel_title <- paste0("Chr",chr," ",mac, " Relative Mutation Rate")
main_rel_out <- paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop.png")

main_scale_title <- paste0("Chr",chr," ",mac, " Mutations--scaled proportions")
main_scale_out <- paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop2.png")

main_dist_title <- paste0("Chr",chr," ",mac, " Distribution by Mutation Type")
main_dist_out <- paste0(imgdir,"/chr",chr,"_",mac,"_dist_count.png")

main_dist_title2 <- paste0("Distribution of Chr",chr," singletons (low exon density bins)")
main_dist_out2 <- paste0(imgdir,"/chr",chr,"_",mac,"_dist_count_all.png")

count_heat_title <- paste0("Chr",chr," ",mac, " Mutation Counts Heatmap")
count_out <- paste0(imgdir,"/chr",chr,"_",mac,"_count_heatmap.png")

gc_heat_title <- paste0("Chr",chr," ",mac," GC Content Heatmap")
gc_heat_out <- paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_gc_heatmap.png")

rel_heat_title <- paste0("Chr",chr,": Relative Mutation Rate Heatmap")
rel_rate_out <- paste0(imgdir,"/chr",chr,"_",mac,"_rel_rate_heatmap.png")

at_seq_title <- paste0("Chr",chr,": AT Mutations by Local Sequence")
at_seq_out <- paste0(imgdir,"/chr",chr,"_",mac,"_AT_seq.png")

gc_seq_title <- paste0("Chr",chr,": GC Mutations by Local Sequence")
gc_seq_out <- paste0(imgdir,"/chr",chr,"_",mac,"_GC_seq.png")

at_rel_prop_title <- paste0("Chr",chr,": AT Relative Mutation Rate by Local Sequence")
at_rel_out <- paste0(imgdir,"/chr",chr,"_",mac,"_AT_rel.png")

gc_rel_prop_title <- paste0("Chr",chr,": GC Relative Mutation Rate by Local Sequence")
gc_rel_out <- paste0(imgdir,"/chr",chr,"_",mac,"_GC_rel.png")

at_map_out <- paste0(imgdir,"/chr",chr,"_",mac,"_AT_map.png")
gc_map_out <- paste0(imgdir,"/chr",chr,"_",mac,"_GC_map.png")

a_map_out <- paste0(imgdir,"/chr",chr,"_",mac,"_A_map.png")
g_map_out <- paste0(imgdir,"/chr",chr,"_",mac,"_G_map.png")
c_map_out <- paste0(imgdir,"/chr",chr,"_",mac,"_C_map.png")
t_map_out <- paste0(imgdir,"/chr",chr,"_",mac,"_T_map.png")