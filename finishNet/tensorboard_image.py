# WoodPixel - Supplementary code for Computational Parquetry:
#             Fabricated Style Transfer with Wood Pixels
#             ACM Transactions on Graphics 39(2), 2020
#
# Copyright (C) 2020 Julian Iseringhausen <opensource@iseringhausen.graphics>
# Copyright (C) 2020 Matthias Hullin, University of Bonn <hullin@cs.uni-bonn.de>
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

import cv2
import io
import os
import keras
import numpy as np
import tensorflow as tf
from utils import load_image_folder, predict_image, compute_layer_shapes


class TensorBoardImage(keras.callbacks.Callback):
    def __init__(self, data_dir, log_dir, config):
        super().__init__()
        self.log_dir = log_dir

        input_shapes, output_shapes = compute_layer_shapes(config)
        self.patch_size_in = input_shapes[0]
        self.patch_size_out = output_shapes[-1]

        self.inputs, self.filenames_in = load_image_folder(
            os.path.join(data_dir, 'input'))

        outputs, filenames_out = load_image_folder(
            os.path.join(data_dir, 'output'))
        summaries_out = [
            self.to_tensorflow_summary(im, tag) for im, tag in zip(
                outputs, filenames_out)]
        self.write_summaries(summaries_out, 0)

    def on_epoch_end(self, epoch, logs={}):
        summaries = []
        for im, f in zip(self.inputs, self.filenames_in):
            im_pred = predict_image(im, self.model)
            summary = self.to_tensorflow_summary(im_pred, f)
            summaries.append(summary)
        self.write_summaries(summaries, epoch)

    def write_summary(self, summary, epoch=0):
        writer = tf.summary.FileWriter(self.log_dir)
        writer.add_summary(summary, epoch)
        writer.close()

    def write_summaries(self, summaries, epoch=0):
        writer = tf.summary.FileWriter(self.log_dir)
        for summary in summaries:
            writer.add_summary(summary, epoch)
        writer.close()

    def to_tensorflow_summary(self, image, tag):
        image = self.to_tensorflow_image(image)
        return tf.Summary(
            value=[tf.Summary.Value(tag=tag, image=image)])

    def to_tensorflow_image(self, tensor):
        height, width, channel = tensor.shape
        image = 255 * tensor
        image = image.astype(np.uint8, copy=False)
        _, image = cv2.imencode('.png', image)
        image_string = image.tostring()
        return tf.Summary.Image(
            height=height, width=width, colorspace=channel,
            encoded_image_string=image_string)
