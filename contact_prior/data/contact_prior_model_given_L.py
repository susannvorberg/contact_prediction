#!/usr/bin/env python

#Compute a prior for a residue pair contact
#based on the observation that the number of contacts increases with
#protein length (dependent on contact_thr and seq_sep)

import numpy as np


contact_prior_model_givenL = {
    6:{
        4:  lambda x:  -0.173 + 0.195 * np.log(x),
        8:  lambda x:  -0.425 + 0.211 * np.log(x),
        12: lambda x:  -0.558 + 0.225 * np.log(x),
        24: lambda x:  -0.741 + 0.233 * np.log(x)
    },
    8:{
        4:  lambda x:  -0.49  + 0.568 * np.log(x),
        8:  lambda x:  -1.226 + 0.598 * np.log(x),
        12: lambda x:  -1.587 + 0.635 * np.log(x),
        24: lambda x:  -2.102 + 0.66  * np.log(x)
    },
    10:{
        4:  lambda x:  -2.017 + 1.429 * np.log(x),
        8:  lambda x:  -3.336 + 1.454 * np.log(x),
        12: lambda x:  -3.939 + 1.502 * np.log(x),
        24: lambda x:  -5.138 + 1.576 * np.log(x)
    }
}