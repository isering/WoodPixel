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
import glob
import numpy as np
import keras
import os
import sys


def load_normalized_image(filename):
    # Load image file and normalize to [0, 1].
    print(filename)
    image = cv2.imread(filename, -1)
    image_dtype = image.dtype
    image = image.astype(np.float32, copy=False)
    if image_dtype == np.uint8:
        image /= 255.0
    elif image_dtype == np.uint16:
        image /= 65535.0
    return image


def load_normalized_mask(filename):
    # Load mask image and noralize to [0,1].
    print(filename)
    image = cv2.imread(filename, cv2.IMREAD_COLOR).astype(
        np.float32, copy=False)
    image[image > 0.0] = 1.0
    return image


def get_image_filenames(folder):
    # Get filenames with valid image extensions.
    image_extensions = ['*.png', '*.tiff', '*.tif', '*.jpg', '*.jpeg']
    filenames = []
    for ext in image_extensions:
        filenames.extend(glob.glob(os.path.join(folder, ext)))
    return sorted(filenames)


def load_image_folder(folder):
    filenames = get_image_filenames(folder)
    images = [load_normalized_image(f) for f in filenames]
    return images, filenames


def load_mask_folder(folder):
    filenames = get_image_filenames(folder)
    masks = [load_normalized_mask(f) for f in filenames]
    return masks, filenames


def load_data_folder(folder):
    # Load input image folder.
    input_images, _ = load_image_folder(os.path.join(folder, 'input'))
    output_images, _ = load_image_folder(os.path.join(folder, 'output'))
    mask_images, _ = load_mask_folder(os.path.join(folder, 'mask'))

    return input_images, output_images, mask_images


def load_input_data(config):
    # Load input images.
    train_input, train_output, train_mask = load_data_folder(
        config['train_dir'])
    test_input, test_output, test_mask = load_data_folder(config['test_dir'])

    # Mask out background pixels.
    train_input = [im * mask for im, mask in zip(train_input, train_mask)]
    train_output = [im * mask for im, mask in zip(train_output, train_mask)]
    test_input = [im * mask for im, mask in zip(test_input, test_mask)]
    test_output = [im * mask for im, mask in zip(test_output, test_mask)]

    if not train_input or not train_output:
        sys.exit('Error: Could not load training images.')

    if not test_input or not test_output:
        sys.exit('Error: Could not load test images.')

    return train_input, train_output, test_input, test_output


def update_progress(job_title, progress):
    # Print progress bar.
    length = 20  # Modify this to change the length.
    block = int(round(length*progress))
    msg = "\r{0}: [{1}] {2}%".format(
        job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
    if progress >= 1:
        msg += " DONE\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()


def save_image_clipped(filename, image):
    image = np.clip(image, 0.0, 1.0)
    image *= 255.0
    image = image.astype(np.uint8, copy=False)
    cv2.imwrite(filename, image)


def predict_image_patch(image, image_out, model, x, y, input_patch_size,
                        output_patch_size, crop_size):
    input_data = np.expand_dims(
        image[y:y+input_patch_size[0], x:x+input_patch_size[1], :], axis=0)
    output_data = np.squeeze(model.predict(input_data))
    image_out[y+crop_size[0]:y+crop_size[0]+output_patch_size[0],
              x+crop_size[1]:x+crop_size[1]+output_patch_size[1],
              :] = output_data


def predict_image(image, model):
    """
    Piecewise predict image.
    """
    image_size = image.shape
    input_patch_size = np.asarray(model.get_layer('layer_in').input_shape[1:3])
    output_patch_size = np.asarray(
        model.get_layer('layer_out').output_shape[1:3])
    crop_size = (input_patch_size - output_patch_size) // 2
    valid_image_size = image_size[:2] - input_patch_size

    image_out = np.zeros(image_size)
    for y in range(0, valid_image_size[0], output_patch_size[0]):
        for x in range(0, valid_image_size[1], output_patch_size[1]):
            predict_image_patch(
                image, image_out, model, x, y, input_patch_size,
                output_patch_size, crop_size)

    for y in range(0, valid_image_size[0], output_patch_size[0]):
        predict_image_patch(
            image, image_out, model, image_size[1] - input_patch_size[1], y,
            input_patch_size, output_patch_size, crop_size)

    for x in range(0, valid_image_size[1], output_patch_size[1]):
        predict_image_patch(
            image, image_out, model, x, image_size[0] - input_patch_size[0],
            input_patch_size, output_patch_size, crop_size)

    predict_image_patch(
        image, image_out, model, image_size[1] - input_patch_size[1],
        image_size[0] - input_patch_size[0], input_patch_size,
        output_patch_size, crop_size)

    return image_out


def compute_layer_shapes(config):
    """
    Compute input layer shape for model with valid padding.
    """
    depth = config['model_depth']
    conv_per_layer = config['conv_per_layer']
    shape_narrow = (config['shape_narrow'], config['shape_narrow'])
    kernel_size = (config['kernel_size'], config['kernel_size'])

    add_boundary = np.asarray(kernel_size) - 1

    input_shapes = []
    output_shapes = []

    current_shape = np.asarray(shape_narrow)
    for _ in range(depth):
        current_shape = current_shape * 2
        output_shapes.insert(0, current_shape)
        current_shape = current_shape + conv_per_layer * add_boundary
        input_shapes.insert(0, current_shape)

    current_shape = np.asarray(shape_narrow)
    for _ in range(depth+1):
        input_shapes.append(current_shape)
        current_shape = current_shape - conv_per_layer * add_boundary
        output_shapes.append(current_shape)
        current_shape = current_shape * 2

    # for in_shape, out_shape in zip(input_shapes, output_shapes):
    #     print('{} -> {}'.format(in_shape, out_shape))

    return input_shapes, output_shapes
