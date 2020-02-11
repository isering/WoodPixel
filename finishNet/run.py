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

import argparse
import cv2
import numpy as np

from keras.models import load_model

from utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True,
                    help='Input image for prediction')
parser.add_argument('--reference', type=str, default='',
                    help='Reference ground truth for comparison')
parser.add_argument('--mask', type=str, default='',
                    help='Input mask')
parser.add_argument('--output', type=str, required=True,
                    help='Output directory')
parser.add_argument('--checkpoint', type=str, required=True,
                    help='Model checkpoint')
parser.add_argument('--patch_size', type=int, default=128,
                    help='Patch size for inference')
FLAGS = parser.parse_args()


if __name__ == '__main__':
    psize = FLAGS.patch_size

    # Load model checkpoint.
    model = load_model(FLAGS.checkpoint)
    model.summary()

    # Load input image and mask.
    inp = load_normalized_image(FLAGS.input)
    if FLAGS.mask:
        mask = load_normalized_mask(FLAGS.mask)
        inp *= mask

    # Predict output.
    patches_y = inp.shape[0] // psize
    patches_x = inp.shape[1] // psize
    image = np.ndarray((patches_y * psize, patches_x * psize, 3), np.float32)
    image = predict_image(inp, model)

    # Write output to disk.
    if not os.path.exists(FLAGS.output):
        os.makedirs(FLAGS.output, )
    save_image_clipped(os.path.join(FLAGS.output, 'input.png'), inp)
    save_image_clipped(os.path.join(FLAGS.output, 'predicted.png'), image)
    if FLAGS.reference:
        reference = load_normalized_image(FLAGS.reference)
        save_image_clipped(
            os.path.join(FLAGS.output, 'reference.png'), reference)
