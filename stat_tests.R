############# TEST PROPORTIONS--GC
print("Test proportions--GC")
aggseq_cg <- aggseq_g[grep("G$", aggseq_g$Category),]
print(prop.test(aggseq_cg$freq, aggseq_cg$COUNT))

aggseq_ta <- aggseq_g[grep("A$", aggseq_g$Category),]
print(prop.test(aggseq_ta$freq, aggseq_ta$COUNT))

aggseq_at <- aggseq_g[grep("T$", aggseq_g$Category),]
print(prop.test(aggseq_at$freq, aggseq_at$COUNT))

############# TEST PROPORTIONS--AT
print("Test proportions--AT")
aggseq_cg <- aggseq_a[grep("G$", aggseq_a$Category),]
print(prop.test(aggseq_cg$freq, aggseq_cg$COUNT))

aggseq_gc <- aggseq_a[grep("C$", aggseq_a$Category),]
print(prop.test(aggseq_gc$freq, aggseq_gc$COUNT))

aggseq_ta <- aggseq_a[grep("A$", aggseq_a$Category),]
print(prop.test(aggseq_ta$freq, aggseq_ta$COUNT))