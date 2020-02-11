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

from __future__ import print_function

import argparse
import os
import random
import sys
import yaml

from keras import backend as K
from keras.callbacks import ModelCheckpoint, TensorBoard
from keras.layers import Add, AveragePooling2D, Concatenate, Conv2D, Cropping2D
from keras.layers import Input, UpSampling2D
from keras.models import Model, Sequential, load_model
from keras.optimizers import Adam

from generator import *
from utils import *
from tensorboard_image import TensorBoardImage

# Parse arguments.
parser = argparse.ArgumentParser()
parser.add_argument('--checkpoint', type=str, default='',
                    help='Load previously saved checkpoint')
parser.add_argument('--config', type=str, default='config.yaml',
                    help='Config file')
parser.add_argument(
    '--data_tensorboard', type=str, default='',
    help='Directory containing images to be evaluated for TensorBoard')
parser.add_argument('--log', type=str, default='', help='Tensorboard log dir.')
parser.add_argument('--model', type=str, required=True,
                    help='Output path to save model')
parser.add_argument('--seed', type=int, default=0,
                    help='Random number generator seed')
parser.add_argument('--workers', type=int, default=8,
                    help='Number of workers for data generation')
FLAGS = parser.parse_args()


def generate_debug_model():
    model = Sequential()
    model.add(Conv2D(3, (3, 3), input_shape=(None, None, 3),
                     padding='same', activation='relu'))
    return model


def generate_model(input_shape, kernel_size, config):
    depth = config['model_depth']
    num_filters = config['num_filters']
    conv_per_layer = config['conv_per_layer']

    layer_in = Input(shape=(None, None, 3))
    current_layer = layer_in

    layers = []

    for _ in range(depth):
        for _ in range(conv_per_layer):
            current_layer = Conv2D(
                num_filters, kernel_size, strides=(1, 1),
                padding='same', activation='relu')(current_layer)
        layers.append(current_layer)
        current_layer = Conv2D(
            num_filters, kernel_size, strides=(2, 2),
            padding='same', activation='relu')(current_layer)
        num_filters *= 2

    current_layer = Conv2D(num_filters, kernel_size, strides=(1, 1),
                           padding='same', activation='relu')(current_layer)

    for layer_skip_connection in reversed(layers):
        num_filters = num_filters // 2
        current_layer = UpSampling2D(size=(2, 2),
                                     interpolation='nearest')(current_layer)

        if config['residual_network']:
            current_layer = Concatenate()(
                [current_layer, layer_skip_connection])

        for _ in range(conv_per_layer):
            current_layer = Conv2D(
                num_filters, kernel_size, strides=(1, 1),
                padding='same', activation='relu')(current_layer)

    if config['residual_network']:
        current_layer = Concatenate()([layer_in, current_layer])

    layer_out = Conv2D(
        3, kernel_size, strides=(1, 1), padding='same',
        activation=config['final_activation'])(current_layer)

    if config['train_color_residual']:
        layer_out = Add()([layer_in, layer_out])

    model = Model(inputs=layer_in, outputs=layer_out)

    return model


def generate_model_unet(shape_narrow, kernel_size, config):
    """
    Generate UNet model with valid padding.
    """
    depth = config['model_depth']
    num_filters = config['num_filters']
    conv_per_layer = config['conv_per_layer']

    kernel_size = np.asarray(kernel_size)
    add_boundary = kernel_size - 1

    input_shapes, output_shapes = compute_layer_shapes(config)
    current_shape = input_shapes[0]

    layer_in = Input(shape=(input_shapes[0][0], input_shapes[0][1], 3),
                     name='layer_in')
    current_layer = layer_in

    layers = []
    for _ in range(depth):
        for _ in range(conv_per_layer):
            current_layer = Conv2D(
                num_filters, kernel_size, strides=(1, 1),
                padding='valid', activation='relu')(current_layer)
            current_shape -= add_boundary
        layers.append(current_layer)
        current_layer = AveragePooling2D()(current_layer)
        num_filters *= 2
        current_shape //= 2

    for _ in range(conv_per_layer):
        current_layer = Conv2D(
            num_filters, kernel_size, padding='valid',
            activation='relu')(current_layer)

    layers = reversed(layers)

    for d, skip_layer in enumerate(layers):
        current_layer = UpSampling2D()(current_layer)
        if config['residual_network']:
            shape_out = output_shapes[depth-d-1]
            shape_in = input_shapes[depth+d+1]
            crop = (shape_out - shape_in) // 2
            layer_crop = Cropping2D(
                cropping=((crop[0], crop[0]), (crop[1], crop[1])))(skip_layer)
            current_layer = Concatenate()([layer_crop, current_layer])
        num_filters //= 2
        current_shape *= 2
        for _ in range(conv_per_layer):
            current_layer = Conv2D(
                num_filters, kernel_size, padding='valid',
                activation='relu')(current_layer)

    layer_out = Conv2D(
        3, (1, 1), padding='valid', activation=config['final_activation'],
        name='layer_out')(current_layer)

    model = Model(inputs=layer_in, outputs=layer_out)
    model.summary()

    return model


def train(train_input, train_output, test_input, test_output, config):
    shape_narrow = (config['shape_narrow'], config['shape_narrow'])
    kernel_size = (config['kernel_size'], config['kernel_size'])

    input_shapes, output_shapes = compute_layer_shapes(config)

    # Create data generators.
    training_generator = Generator(
        'train', train_input, train_output, config['samples_train'],
        input_shapes[0], output_shapes[-1], config['batch_size'])
    test_generator = Generator(
        'test', test_input, test_output, config['samples_test'],
        input_shapes[0], output_shapes[-1], config['batch_size'],
        random_samples=False, shuffle=False)

    if FLAGS.checkpoint:
        model = load_model(FLAGS.checkpoint)
    else:
        model = generate_model_unet(shape_narrow=shape_narrow,
                                    kernel_size=kernel_size, config=config)
        optimizer = Adam(lr=float(config['learning_rate']))
        model.compile(optimizer=optimizer, loss=config['loss'])

    # Add callbacks.
    callback_list = []

    if not os.path.exists(FLAGS.model):
        os.makedirs(FLAGS.model)
    filepath = os.path.join(FLAGS.model, '{epoch:04d}.hdf5')
    checkpoint = ModelCheckpoint(
        filepath, monitor='val_loss',
        verbose=1, save_best_only=False, mode='auto',
        period=config['save_period'])
    callback_list.append(checkpoint)

    if FLAGS.log:
        tensor_board_callback = TensorBoard(
            log_dir=FLAGS.log, batch_size=config['batch_size'])
        callback_list.append(tensor_board_callback)
        if FLAGS.data_tensorboard:
            tensor_board_image = TensorBoardImage(
                FLAGS.data_tensorboard, FLAGS.log, config)
            callback_list.append(tensor_board_image)

    # Train model on dataset.
    model.fit_generator(generator=training_generator,
                        validation_data=test_generator,
                        epochs=config['epochs'],
                        callbacks=callback_list,
                        use_multiprocessing=False,
                        max_queue_size=50,
                        workers=FLAGS.workers)


def load_config(filename):
    with open(filename, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as e:
            sys.exit(e)
    return config


if __name__ == '__main__':
    # Initialize random number generator.
    random.seed(FLAGS.seed)

    # Load config file.
    config = load_config(FLAGS.config)
    config['train_dir'] = os.path.join(
        os.path.dirname(FLAGS.config), config['train_dir'])
    config['test_dir'] = os.path.join(
        os.path.dirname(FLAGS.config), config['test_dir'])

    # Load input data.
    train_input, train_output, test_input, test_output = load_input_data(
        config)

    # Start training.
    train(train_input, train_output, test_input, test_output, config)
