#ifndef WEIGHTING_H
#define WEIGHTING_H

#include <stdint.h>
#define GAP 20

void count_ids(
	const uint8_t *msa,
	uint64_t *ids,
	const uint64_t nrow,
	const uint64_t ncol
);

void calculate_weights_simple(
	const uint8_t *msa,
	double *weights,
	double cutoff,
	bool ignore_gaps,
	const uint64_t nrow,
	const uint64_t ncol
);


#endif
