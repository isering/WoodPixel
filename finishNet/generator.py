# WoodPixel - Supplementary code for Computational Parquetry:
#             Fabricated Style Transfer with Wood Pixels
#             ACM Transactions on Graphics 39(2), 2020
#
# Copyright (C) 2020  Julian Iseringhausen, University of Bonn, <iseringhausen@cs.uni-bonn.de>
# Copyright (C) 2020  Matthias Hullin, University of Bonn, <hullin@cs.uni-bonn.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This file contains the data generator class.
import keras
import numpy as np
import random

from utils import *

# FIXME(isering): random_samples=True does not work for true multiprocessing,
#                 since the initial state of the Python interpreter (including
#                 the state of the random number generator) is copied to each
#                 process. Thus each process would produce the same samples.
#                 For multi-threading based parallelization this should work
#                 fine, since the threads share the same random number
#                 generator state.


class Generator(keras.utils.Sequence):
    # Generates data for Keras.

    def __init__(self, name, inputs, outputs, num_samples, shape_in, shape_out,
                 batch_size, random_samples=True, shuffle=True,
                 n_channels=3):
        # Initialization.
        self.name = name
        self.inputs = inputs
        self.outputs = outputs
        self.shape_in = shape_in
        self.shape_out = shape_out
        self.crop = (shape_in - shape_out) // 2
        self.batch_size = batch_size
        self.num_samples = num_samples
        self.n_channels = n_channels
        self.random_samples = random_samples
        self.shuffle = shuffle
        if not self.random_samples:
            self.samples = self.generate_samples(self.num_samples)

    def __len__(self):
        # Denotes the number of batches per epoch.
        return self.num_samples // self.batch_size

    def __getitem__(self, index):
        # Get samples.
        if self.random_samples:
            samples_batch = self.generate_samples(self.batch_size)
        else:
            samples_batch = self.samples[
                index*self.batch_size:(index+1)*self.batch_size]

        # Generate data for samples.
        X, Y = self.__data_generation(samples_batch)

        return X, Y

    def generate_samples(self, num_samples):
        # Generate samples.
        num_inputs = len(self.inputs)

        samples = []
        for _ in range(num_samples):
            sample = {}
            sample['image_id'] = random.randrange(num_inputs)
            sample['patch_pos_x'] = random.random()
            sample['patch_pos_y'] = random.random()
            sample['rotation'] = random.randrange(4)
            sample['flip_lr'] = random.randint(0, 1)
            sample['scale'] = random.uniform(0.9, 1.1)
            samples.append(sample)
        return samples

    def on_epoch_end(self):
        # Updates indices after each epoch.
        if not self.random_samples and self.shuffle:
            np.random.shuffle(self.samples)

    def __data_generation(self, samples):
        # Generates data containing batch_size samples.
        X = np.empty((self.batch_size, *self.shape_in, self.n_channels))
        Y = np.empty((self.batch_size, *self.shape_out, self.n_channels))

        for i, sample in enumerate(samples):
            image_shape = self.inputs[sample['image_id']].shape
            pos_y = int(
                (image_shape[0]-self.shape_in[0]-1) * sample['patch_pos_y'])
            pos_x = int(
                (image_shape[1]-self.shape_in[1]-1) * sample['patch_pos_x'])

            x = self.inputs[sample['image_id']][
                pos_y:pos_y+self.shape_in[0], pos_x:pos_x+self.shape_in[1], :]
            y = self.outputs[sample['image_id']][
                pos_y+self.crop[0]:pos_y+self.crop[0]+self.shape_out[0],
                pos_x+self.crop[1]:pos_x+self.crop[1]+self.shape_out[1], :]

            x = np.rot90(x, sample['rotation'])
            y = np.rot90(y, sample['rotation'])

            if sample['flip_lr']:
                x = np.fliplr(x)
                y = np.fliplr(y)

            x = x * sample['scale']
            y = y * sample['scale']

            np.clip(x, 0.0, 1.0, out=x)
            np.clip(y, 0.0, 1.0, out=y)

            # cv2.imshow('x', x)
            # cv2.imshow('y', y)
            # cv2.waitKey(0)

            X[i, ] = x
            Y[i, ] = y

        return X, Y
