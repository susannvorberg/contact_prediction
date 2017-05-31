#!/usr/bin/env python

import os
import json
import utils.io_utils as io



class AlignmentFeatures():
    """
    Compute sequence and alignment derived features
    """

    def __init__(self, alignment_file):

        self.alignment_file = alignment_file
        self.protein=os.path.basename(self.alignment_file)
        self.msa = io.read_alignment(alignment_file)

        self.features = {}

    def __repr__(self):
        repr_str ="Features for protein {0}: \n".format(self.protein)
        for key, value in self.features.iteritems():
            repr_str += "{0} : {1}\n".format(key, len(value))

    def read_features(self, feature_file):
        self.features = json.load(feature_file)

    def write_features(self, feature_file):
        json.dump(self.features, feature_file)

    def get_features(self):
        return self.features



