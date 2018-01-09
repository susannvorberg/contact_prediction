#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

#include "weighting.h"

/**
 * Count the number of sequence identities for all rows in an MSA
 *
 * @param[in] seq The MSA to work on
 * @param[out] counts The number of sequence identities
 * @param[in] nrow The number of columns in the MSA
 * @param[in] ncol The number of rows in the MSA
 */
void count_ids(
	const uint8_t *msa,
	uint64_t *ids,
	const uint64_t nrow,
	const uint64_t ncol
) {
	uint64_t nij = nrow * (nrow + 1) / 2;

	omp_set_dynamic(0);

	#pragma omp parallel
	{
		uint64_t ij;

		#pragma omp for nowait private(ij)
		for(ij = 0; ij < nij; ij++) {

			// compute i and j from ij
			// http://stackoverflow.com/a/244550/1181102
			uint64_t i, j;
			{
				uint64_t ii = nrow * (nrow + 1) / 2 - 1 - ij;
				uint64_t K = floor((sqrt(8 * ii + 1) - 1) / 2);
				i = nrow - 1 - K;
				j = ij - nrow * i + i * (i + 1) / 2;
			}

			uint64_t my_ids = 0;
			for(uint64_t k = 0; k < ncol; k++) {
				if(msa[i * ncol + k] == msa[j * ncol + k]) {
					my_ids++;
				}
			}

			ids[i * nrow + j] = my_ids;
		}
	}
}


void calculate_weights_simple(
	const uint8_t *msa,
	double *weights,
	double cutoff,
	bool ignore_gaps,
	const uint64_t nrow,
	const uint64_t ncol
) {
	uint64_t nij = nrow * (nrow + 1) / 2;

	omp_set_dynamic(0);

	#pragma omp parallel
	{

		uint64_t ij;

		#pragma omp for nowait private(ij)
		for(ij = 0; ij < nij; ij++) {

			// compute i and j from ij
			// http://stackoverflow.com/a/244550/1181102
			uint64_t i, j;
			{
				uint64_t ii = nrow * (nrow + 1) / 2 - 1 - ij;
				uint64_t K = floor((sqrt(8 * ii + 1) - 1) / 2);
				i = nrow - 1 - K;
				j = ij - nrow * i + i * (i + 1) / 2;
			}


			uint64_t my_ids = 0;
			uint64_t idthres = ceil(cutoff * ncol);

			if (ignore_gaps){
				uint64_t ncol_ij = ncol;
				for(uint64_t k = 0; k < ncol; k++) {
					if(msa[i * ncol + k] == msa[j * ncol + k] ) {
						if(msa[i * ncol + k] == GAP) ncol_ij--;
						else my_ids++;
					}
				}
				idthres = ceil(cutoff * ncol_ij);

			}else{
				for(uint64_t k = 0; k < ncol; k++) {
					if(msa[i * ncol + k] == msa[j * ncol + k] ) {
						my_ids++;
					}
				}

			}


			if(my_ids >= idthres) {
				#pragma omp atomic
				weights[i]++;
				#pragma omp atomic
				weights[j]++;
			}

		}
	}

	for(uint64_t i = 0; i < nrow; i++) {
		weights[i] = 1.0 / (weights[i] - 1);
	}

	fflush(stdout);
}
