kmc_histogram: kmc_histogram.cpp
	g++ -std=c++17 -O3 -I../external/KMC/kmc_api \
	kmc_histogram.cpp \
	../external/KMC/kmc_api/kmc_file.cpp \
	../external/KMC/kmc_api/kmc_file.h \
	../external/KMC/kmc_api/mmer.cpp \
	../external/KMC/kmc_api/mmer.h \
	../external/KMC/kmc_api/kmer_api.cpp \
	../external/KMC/kmc_api/kmer_api.h \
	../external/KMC/kmc_api/kmer_defs.h \
	-lpthread -o kmc_histogram

clean:
	rm -f kmc_histogram