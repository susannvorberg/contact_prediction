#!/usr/bin/env python

from keras.models import Sequential
from keras.layers import Dense




alignment_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/psicov/"
pdb_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/pdb_renum_combs/"
psipred_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/psipred/hhfilter_results_n5e01/"
netsurfp_dir="/home/vorberg/work/data/benchmarkset_cathV4.1/netsurfp/"
property_files_dir = "/home/vorberg/work/data/benchmarkset_cathV4.1/dataset/dataset_properties/"
nr_contacts=1000
nr_non_contacts=5000
braw_dir=None
window_size = 5
seq_separation          = 8
contact_threshold       = 8
non_contact_threshold   = 25




# create model
model = Sequential()
model.add(Dense(feature_df.shape[1], input_dim=feature_df.shape[1], kernel_initializer='normal', activation='relu'))
model.add(Dense(1, kernel_initializer='normal', activation='sigmoid'))
# Compile model, logarithmic loss function = binary_crossentropy
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])




# x_train and y_train are Numpy arrays --just like in the Scikit-Learn API.
model.fit(feature_df.as_matrix(), class_df['contact'].values, epochs=20, batch_size=128)


score = model.evaluate(feature_df.as_matrix(), class_df['contact'].values, batch_size=128)
print("\n%s: %.2f%%" % (model.metrics_names[1], score[1]*100))