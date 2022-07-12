import os
import sys
import copy
import scipy.stats
import argparse
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import layers
from matplotlib import pyplot as plt

# if tf.test.gpu_device_name():
#     print('Default GPU Device:{}'.format(tf.test.gpu_device_name()))
#     print("It should print Default GPU Device: /device:GPU:0")
# else:
#     print("Please install GPU version of TF")

print(tf.version.VERSION)
print(sys.version)

# physical_devices = tf.config.experimental.list_physical_devices('GPU')
# if len(physical_devices) > 0:
#     tf.config.experimental.set_memory_growth(physical_devices[0], True)

def get_train_and_test_splits(dataset, train_size, settings=argparse.Namespace()):
    # Load relevant data from data file and parse it into paired (N, k) and (N, ) feature and score matrices
    lines = open(dataset, 'r').readlines()[1:]
    history = argparse.Namespace()
    history.muts = [line.split(':')[0].split(',') for line in lines]
    history.score = [float(line.split(':')[1]) for line in lines]

    # Gather thread history.muts and history.score attributes from each thread together into respective lists
    muts = []
    scores = []
    for ii in range(len(history.score)):
        if history.score[ii]:
            scores.append(history.score[ii])
            muts.append(history.muts[ii])     # only get muts with corresponding scores
        else:
            print(history.muts[ii])
    print('number of distinct mutants: ' + str(len(muts)))
    # print(muts)

    # Normalize scores to between 0 and 1
    print('min: ' + str(min(scores)))
    print('max: ' + str(max(scores)))
    scores = [(score - min(scores)) / (max(scores) - min(scores)) for score in scores]

    # Convert each mutants list entry into full sequence-length list of integer-encoded features
    feats = []  # initialize feature vector
    if settings.encoding == 'one_hot' or settings.encoding == 'integer':
        for mut in muts:
            feats.append(integer_encoder(mut, settings.seq))    # get integer encoding
        if settings.encoding == 'one_hot':
            feats = keras.utils.to_categorical(feats)           # convert integer encoding to one-hot encoding
    elif settings.encoding == 'embedded':
        feats = embedded_encoder(muts, settings.seq)  # learn an embedding
    elif settings.encoding == 'categorical':
        for mut in muts:
            feats.append(categorical_encoder(mut, settings.seq))
    else:
        raise RuntimeError('unrecognized encoding: ' + settings.encoding)

    # Subdivide into training and testing datasets based on train_size argument and return
    rng = np.random.default_rng()
    indices = [ind for ind in range(len(feats))]
    print("indices length:" + str(len(indices)))
    rng.shuffle(indices)
    test_split = (len(feats) - train_size) / len(feats)


    ## A kludge for now that allows me to specify mutants to withhold for the test set
    # pre_indices = copy.copy(indices)
    # to_withold = [['64ASP', '386SER', '394ALA'],['64ASP', '386SER'],['64ASP', '394ALA', '407PRO', '423PRO'],['64ASP', '386SER', '360VAL', '389ASN']]#,['64ASP', '386SER', '422THR']]
    # indices_of_witholding = [muts.index(mut) for mut in to_withold] # order of single muts within each mut matters for this implementation
    # print(indices_of_witholding)
    # for index in indices_of_witholding:
    #     print(muts[index])
    # to_exchange = []
    # for ii in indices_of_witholding:
    #     if indices.index(ii) < int(len(feats) * (1 - test_split)):  # if ii is in the train split
    #         to_exchange.append(ii)
    # for jj in range(len(to_exchange)):
    #     temp = indices[-1 * (jj + 1)]
    #     original_index = indices.index(to_exchange[jj])
    #     indices[-1 * (jj + 1)] = to_exchange[jj]
    #     indices[original_index] = temp
    #     print('moving mutant at index ' + str(to_exchange[jj]) + ' (' + str(muts[to_exchange[jj]]) + ') to index: ' + str(-1 * (jj + 1)))

    muts = [muts[ii] for ii in indices]
    feats = [feats[ind] for ind in indices]
    scores = [scores[ind] for ind in indices]

    # assert (all([item in indices for item in pre_indices]))
    # assert (all([mut in muts[int(len(feats) * (1 - test_split)):] for mut in to_withold]))

    feats_train = np.array(feats[:int(len(feats) * (1 - test_split))])
    scores_train = np.array(scores[:int(len(feats) * (1 - test_split))])
    feats_test = np.array(feats[int(len(feats) * (1 - test_split)):])
    scores_test = np.array(scores[int(len(feats) * (1 - test_split)):])
    #assert all([keras.utils.to_categorical(integer_encoder(mut, settings.seq)) in feats_test for mut in to_withold])

    train_dataset = tf.data.Dataset.from_tensor_slices((feats_train, scores_train))
    test_dataset = tf.data.Dataset.from_tensor_slices((feats_test, scores_test))

    return train_dataset, test_dataset

def integer_encoder(muts, wt_seq):
    # Encode sequence with mutations as a list of integers

    all_resnames = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA',
                    'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']    # all resnames in encoding order

    if muts == ['WT']:   # special case
        muts = []        # no mutation applied

    values = copy.copy(wt_seq)
    # print(wt_seq)
    # print(values)
    # print(muts)
    # print("muts")
    # print(muts)
    # print(type(muts))
    index = 0
    for mut in muts:
        print(mut)
        # print(int(mut[:-3]) - 1)
        # print(type(mut))
        # index = index + 1
        # print(index)
        # print(mut)
        values[int(mut[:-3]) - 1] = mut[-3:]


    int_encoded = []
    values_index = -1
    for item in values:

        int_encoded += [all_resnames.index(item)]

        # values_index += 1
        # if wt_seq[values_index] == item:
        #     int_encoded += [0]
        # else:
        #     int_encoded += [all_resnames.index(item) + 1]

    return int_encoded

def categorical_encoder(muts, wt_seq):
    # Encode sequence with mutations as a list of resnames

    all_resnames = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA',
                    'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']    # all resnames in enconding order

    if muts == ['WT']:   # special case
        muts = []        # no mutation applied

    values = copy.copy(wt_seq)
    for mut in muts:
        values[int(mut[:-3]) - 1] = mut[-3:]

    return values

def embedded_encoder(muts, wt_seq):
    # Learn an embedding for the given mutants in terms of the following features:
    #   Residue index
    #   Geometric positions (x, y, z coordinates) of alpha carbons  # todo: implement
    #   Sidechain charge at pH 7
    #   Molecular weight
    #   Sidechain hydrophobicity
    #   Sidechain aromaticity

    embedding_size = 1

    # First, define by-residue-name characteristics
    resnames = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA',
                'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']    # all resnames in enconding order
    charges =  [1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]   # charges in same order
    mws =      [174, 155, 146, 133, 146, 105, 119, 132, 146, 121, 75, 115, 89, 117, 131, 131, 149, 165, 181, 204]
    hydrophobicity = [0, 0.165, 0.283, 0.028, 0.043, 0.359, 0.450, 0.236, 0.251, 0.680, 0.501, 0.711, 0.616,
                      0.825, 0.943, 0.943, 0.738, 1, 0.880, 0.878]
    aromaticity = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]

    def vocabularize(raw):
        # Convert a "raw" list of values into a categorical vocabulary for use in Embedding
        vocab_index = 0
        vocab = []
        for ii in range(len(raw)):
            if raw[ii] in raw[:ii]:
                vocab.append(vocab[raw.index(raw[ii])])
            else:
                vocab.append(vocab_index)
                vocab_index += 1
        return vocab

    # Build an Input for each feature
    charge_input = keras.Input(name='charge', shape=[441])
    mw_input = keras.Input(name='mw', shape=[441])
    hydrophobicity_input = keras.Input(name='hydrophobicity', shape=[1])
    aromaticity_input = keras.Input(name='aromaticity', shape=[1])
    # coords_input = keras.Input(name='coords', shape=[3])

    # Build an Embedding layer for each Input
    charge_embedding = keras.layers.Embedding(name='charge_embedding',
                               input_dim=max(vocabularize(charges)) + 1,
                               input_length=len(wt_seq),
                               output_dim=embedding_size)(charge_input)
    mw_embedding = keras.layers.Embedding(name='mw_embedding',
                               input_dim=max(vocabularize(mws)) + 1,
                               input_length=len(wt_seq),
                               output_dim=embedding_size)(mw_input)
    # hydrophobicity_embedding = keras.layers.Embedding(name='hydrophobicity_embedding',
    #                            input_length=len(wt_seq),
    #                            output_dim=embedding_size)(hydrophobicity_input)
    # aromaticity_embedding = keras.layers.Embedding(name='aromaticity_embedding',
    #                            input_length=len(wt_seq),
    #                            output_dim=embedding_size)(aromaticity_input)
    # coords_embedding = keras.layers.Embedding(name='coords_embedding',
    #                            input_dim=len(resnames),
    #                            output_dim=embedding_size)(coords_input)

    # Merge the layers with a dot product along the second axis
    merged = keras.layers.Dot(name='dot_product', normalize=True, axes=2)([charge_embedding, mw_embedding])#, hydrophobicity_embedding, aromaticity_embedding])

    # Reshape to be a single number (shape will be (None, 1))
    # merged = keras.layers.Reshape(target_shape=[1])(merged)

    # Output neuron
    out = keras.layers.Dense(1, activation='sigmoid')(merged)
    model = keras.Model(inputs=[charge_input, mw_input], outputs=out)
    # model = keras.Model(inputs=[charge_input, mw_input, hydrophobicity_input, aromaticity_input], outputs=out)

    # Minimize binary cross entropy
    model.compile(optimizer='Adam', loss='binary_crossentropy', metrics=['accuracy'])

    print(model.summary())

    charges = vocabularize(charges)
    mws = vocabularize(mws)

    # Featurize each item in muts
    featurized_muts = []
    for mut in muts:
        # Build full mutated sequence for this mut
        values = copy.copy(wt_seq)
        for mu in mut:
            values[int(mu[:-3]) - 1] = mu[-3:]

        featurized_muts.append([[charges[resnames.index(res)] for res in values],
                                [mws[resnames.index(res)] for res in values]])

    featurized_muts = np.array(featurized_muts)
    featurized_muts = {'charge': featurized_muts[:, 0], 'mw': featurized_muts[:, 1]}

    embeddings = model.predict(featurized_muts)

    return embeddings


"""
## Compile, train, and evaluate the model
"""

def run_experiment(model, loss, train_dataset, test_dataset, num_epochs):

    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=0.0001),
        # optimizer='adam',
        loss=loss,
        metrics=[keras.metrics.MeanSquaredError()],
        # metrics=tf.keras.metrics.categorical_accuracy
    )

    # callbacks = [
    #     EarlyStopping(patience=15, min_delta=0.00001)
    # ]

    print("Start training the model...")
    history = model.fit(train_dataset, epochs=num_epochs, validation_data=test_dataset)
    print("Model training finished.")
    _, rmse = model.evaluate(train_dataset, verbose=0)
    print(f"Train RMSE: {round(rmse, 3)}")

    print("Evaluating model performance...")
    _, rmse = model.evaluate(test_dataset, verbose=0)
    print(f"Test RMSE: {round(rmse, 3)}")

    return history

def create_model(batch_size, length):
    input_shape = (batch_size, length, 20)     # set appropriate input size
    print("batch size")
    print(batch_size)
    print("length")
    print(length)

    # pab1 model
    # model = tf.keras.Sequential([
    #        tf.keras.layers.Input(batch_input_shape=input_shape),
    # #       # tf.keras.layers.InputLayer(batch_input_shape=input_shape),
    # #       # tf.expand_dims(batch_size, axis=1),
    # #       # tf.keras.layers.Embedding(20, 15),
    # #       # tf.expand_dims(axis=-1),
    #        tf.keras.layers.Conv1D(128, 17, strides=1, padding='valid', activation=tf.nn.leaky_relu, input_shape=input_shape),  # relu performs much better than linear
    #        tf.keras.layers.Conv1D(128, 17, strides=1, padding='valid', activation=tf.nn.leaky_relu),
    # #       # tf.keras.layers.MaxPool1D(2),
    # #
    # #       # tf.keras.layers.BatchNormalization(1),      # subjectively performing much better than dropout
    #        tf.keras.layers.Conv1D(128, 17, strides=1, padding='valid', activation=tf.nn.leaky_relu),
    # #       # tf.keras.layers.MaxPool1D(2),
    # #
    # #       # tf.keras.layers.BatchNormalization(1),
    # #       tf.keras.layers.Conv1D(128, 17, strides=1, padding='valid', activation=tf.nn.leaky_relu),
    # #       # tf.keras.layers.MaxPool1D(2),
    # #       # tf.keras.layers.GlobalAveragePooling1D(),   # global avg. pooling or flatten --> dense "head" network
    #        tf.keras.layers.Flatten(),
    # #
    #        tf.keras.layers.Dense(100, activation=tf.nn.leaky_relu),
    #        tf.keras.layers.Dropout(0.2),
    #        tf.keras.layers.Dense(1)
    # ])

    # bgl3 model
    model = tf.keras.Sequential([
          # tf.keras.layers.Input(batch_input_shape=input_shape),
          tf.keras.layers.Input(input_shape=input_shape),
          tf.keras.layers.Conv1D(128, 3, strides=1, padding='valid', activation=tf.nn.leaky_relu),
          tf.keras.layers.Conv1D(128, 3, strides=1, padding='valid', activation=tf.nn.leaky_relu),
          tf.keras.layers.Conv1D(128, 3, strides=1, padding='valid', activation=tf.nn.leaky_relu),

          tf.keras.layers.Flatten(),

          tf.keras.layers.Dense(100, activation=tf.nn.leaky_relu),
          tf.keras.layers.Dropout(0.2),
          tf.keras.layers.Dense(1)
    ])

    # ube4b model
    # model = tf.keras.Sequential([
    #       tf.keras.layers.Input(batch_input_shape=input_shape),
    #       tf.keras.layers.Conv1D(32, 17, strides=1, padding='valid', activation=tf.nn.leaky_relu),
    #
    #       tf.keras.layers.Flatten(),
    #
    #       tf.keras.layers.Dense(100, activation=tf.nn.leaky_relu),
    #       tf.keras.layers.Dropout(0.2),
    #       tf.keras.layers.Dense(1)
    #   ])

    model.summary()

    return model

def main(dataset, settings):
    best_r = None
    for training_index in range(max(1, settings.repeat_trainings)):
        full_size = len(open(dataset, 'r').readlines()) - 1
        print("full_size:" + str(full_size))
        train_size = int(0.80 * full_size)

        train_dataset, test_dataset = get_train_and_test_splits(dataset, train_size, settings)

        num_epochs = 100    # should choose a number such that the training loss has just leveled out at the end
        BATCH_SIZE = 32     # 64 ube4b, 128 pab1/bgl3 bgl3 32 for test
        # Implement random_forest support (ignore if you're not using random forest)
        if settings.model_type == 'random_forest':
            import tensorflow_decision_forests as tfdf
            from wurlitzer import sys_pipes
            import pandas as pd
            model = tfdf.keras.RandomForestModel(task=tfdf.keras.Task.REGRESSION, num_trees=300)
            model.compile(metrics=["mse"])

            # Convert datasets to pandas dataframe
            ii = 0
            with open('train.csv', 'w') as f:
                f.write('score,' + ','.join([str(ii) for ii in range(441)]))
                f.write('\n')
                for item in train_dataset:
                    f.write(str(tf.keras.backend.get_value(item[1])) + ',' + ','.join([item.decode() for item in tf.keras.backend.get_value(item[0])]))
                    f.write('\n')
                    ii += 1
            train_dataset = pd.read_csv('train.csv', header=0)
            train_dataset = tfdf.keras.pd_dataframe_to_tf_dataset(train_dataset, label='score', task=tfdf.keras.Task.REGRESSION)
            with open('test.csv', 'w') as f:
                f.write('score,' + ','.join([str(ii) for ii in range(441)]))
                f.write('\n')
                for item in test_dataset:
                    f.write(str(tf.keras.backend.get_value(item[1])) + ',' + ','.join([item.decode() for item in tf.keras.backend.get_value(item[0])]))
                    f.write('\n')
                    ii += 1
            test_dataset = pd.read_csv('test.csv', header=0)
            test_dataset = tfdf.keras.pd_dataframe_to_tf_dataset(test_dataset, label='score', task=tfdf.keras.Task.REGRESSION)

            model.fit(train_dataset)
            # os.remove('temp.csv')
            # with sys_pipes():
            #     model.fit(x=train_dataset)
            model.summary()
            evaluation = model.evaluate(test_dataset, return_dict=True)
            print(evaluation)
            print(f"MSE: {evaluation['mse']}")

            # Get samples to use in validation
            sample = full_size - train_size  # number of samples
            examples, targets = list(test_dataset.unbatch().shuffle(BATCH_SIZE * 10).batch(sample))[0]

            predicted = model(examples).numpy()
            r = scipy.stats.pearsonr([predicted[ii][0] for ii in range(sample)], [targets[ii] for ii in range(sample)])
            predicted = model(examples).numpy()
            for idx in range(sample):
                print(f"Predicted: {round(float(predicted[idx][0]), 3)} - Actual: {round(float(targets[idx]), 3)}")
            print('Pearson r (p): ' + str(r))

            sys.exit()

        # If not random_forest, proceed with CNN
        model = create_model(BATCH_SIZE, len(settings.seq))

        SHUFFLE_BUFFER_SIZE = 100
        train_dataset = train_dataset.shuffle(SHUFFLE_BUFFER_SIZE).batch(BATCH_SIZE)
        test_dataset = test_dataset.batch(BATCH_SIZE)

        # Reload previous model if desired
        if os.path.exists('last_model.keras') and settings.reload_weights:
            model.load_weights('last_model.keras')

        history = run_experiment(model, keras.losses.MeanSquaredError(), train_dataset, test_dataset, num_epochs)

        print("successfully ran experiment!!!")

        # Get samples to use in validation
        sample = full_size - train_size    # number of samples

        print("reached examples and targets line")
        examples, targets = list(test_dataset.unbatch().shuffle(BATCH_SIZE * 10).batch(sample))[0]

        print("successfully assigned examples and targets")

        print("reached predictions part")

        # print(type(model(examples)))

        # predicted = []
        # for ex in examples:
        #     print('model(ex).numpy() type')
        #     print(type(model(ex).numpy()))
        #     predicted.append(model(ex).numpy())


        examples_first_quart = examples[:int(len(examples)/4)]
        examples_second_quart = examples[int(len(examples)/4):int(len(examples)/2)]
        examples_third_quart = examples[int(len(examples)/2):int(len(examples)/(4/3))]
        examples_fourth_quart = examples[int(len(examples)/(4/3)):]

        print("reaches first_pred")
        first_pred = model(examples_first_quart).numpy()

        print("reaches second pred")
        second_pred = model(examples_second_quart).numpy()

        print("reaches third pred")
        third_pred = model(examples_third_quart).numpy()

        print("reaches fourth pred")
        fourth_pred = model(examples_fourth_quart).numpy()

        print(first_pred)
        print(second_pred)
        print(third_pred)
        print(fourth_pred)

        predicted = np.concatenate((first_pred, second_pred, third_pred, fourth_pred), axis=0)
        print(predicted)

        # predicted = model(examples)

        # print("passed model line! on numpy line")

        # predicted = predicted.numpy()

        print("successfully assigned predictions!!!")

        print("reached this_r part")

        this_r = scipy.stats.pearsonr([predicted[ii][0] for ii in range(sample)], [targets[ii] for ii in range(sample)])

        print("successfully passed this_r part!!!")

        if best_r == None or this_r > best_r:
            best_r = this_r
            best_model = model
            best_sample = sample
            best_examples = examples
            best_targets = targets
            best_history = history

    print('Pearson r (p): ' + str(best_r))

    best_model.save('last_model.keras')  # dump model for recovery later if desired

    print(best_history.history)
    # plt.semilogy(best_history.history['loss'])
    # plt.semilogy(best_history.history['val_loss'])
    # plt.xlabel('epoch')
    # plt.ylabel('loss')
    # plt.legend(['train', 'test'], loc='upper left')
    # plt.show()
    open('last_loss_curve_bgl3_ss_2880_lim.txt','w').write(str(best_history.history['loss']))
    open('last_val_loss_curve_bgl3_ss_2880_lim.txt', 'w').write(str(best_history.history['val_loss']))

    predicted = best_model(best_examples).numpy()


    for idx in range(best_sample):
        print(f"Predicted: {round(float(predicted[idx][0]), 3)} - Actual: {round(float(targets[idx]), 3)}")
        # if round(float(targets[idx]), 7) in [round(val, 7) for val in [0.246499653, 0.277258577, 0.292320875, 0.246250577, 0.200327526]]:
        #     colors.append('#ff7f0e')
        # else:
        #     colors.append('#1f77b4')
    print('Pearson r (p): ' + str(best_r))

    pred = [predicted[idx][0] for idx in range(best_sample)]

    open('pred_ube_ss_2880_lim.txt', 'w').write(str(pred))
    open('best_targets_bgl3_ss_2880_lim.txt', 'w').write(str(best_targets))

    # plt.scatter(pred, best_targets, s=3)#, c=colors)
    # plt.plot([0, 1], [0, 1])
    # plt.xlabel('Predicted')
    # plt.ylabel('Actual')
    # plt.show()

if __name__ == "__main__":
    settings = argparse.Namespace()
    settings.reload_weights = False
    settings.encoding = 'one_hot'
    settings.repeat_trainings = 1
    settings.model_type = 'cnn'

    # if tf.test.gpu_device_name():
    #     print('Default GPU Device: {}'.format(tf.test.gpu_device_name()))
    # else:
    #     print("Please install GPU version of TF")
    # print("whether gpu is being used")
    # print(tf.config.list_physical_devices('GPU'))


    # print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

    # os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"

    # gpu_options = tf.GPUOptions(visible_device_list="/device:GPU:0")
    # sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

    # tf.debugging.set_log_device_placement(True)

    # try:
    #     # Specify an invalid GPU device
    #     with tf.device('/device:GPU:0'):
    #         a = tf.constant([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    #         b = tf.constant([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    #         c = tf.matmul(a, b)
    # except RuntimeError as e:
    #     print(e)

    config = tf.compat.v1.ConfigProto()
    config.gpu_options.allow_growth = True
    session = tf.compat.v1.Session(config=config)

    # this stuff is needed
    seq_and_score = input('Enter file name: ')
    # seq_and_score = 'ube4b_MLformat_all_90960.txt'  # set the path to the training data file here
    # seq_and_score = 'ube4b_MLformat_80.txt'
    # seq_and_score = 'bgl_MLformat_all.txt'  # set the path to the training data file here

    # seq_and_score = 'pab1_MLformat_ss_2880_lim.txt' # predictor_clean.py_pab1_ss_2880_lim
    # seq_and_score = 'pab1_MLformat_not_ss_2880_lim.txt' # predictor_clean.py_pab1_not_ss_2880_lim
    # seq_and_score = 'bgl3_MLformat_ss_2880_lim.txt' # predictor_clean.py_bgl3_ss_2880_lim
    # seq_and_score = 'bgl3_MLformat_not_ss_2880_lim.txt' # predictor_clean.py_bgl3_not_ss_2880_lim
    # seq_and_score = 'ube4b_MLformat_ss_2880_lim.txt' # predictor_clean.py_ube4b_ss_2880_lim
    # seq_and_score = 'ube4b_MLformat_not_ss_2880_lim.txt' # predictor_clean.py_ube4b_not_ss_2880_lim

    # seq_and_score = 'pab1_MLformat_ss_v2_2880_lim.txt' # predictor_clean.py_pab1_ss_v2_2880_lim
    # seq_and_score = 'pab1_MLformat_not_ss_v2_2880_lim.txt' # predictor_clean.py_pab1_not_ss_v2_2880_lim
    # seq_and_score = 'bgl3_MLformat_ss_v2_2880_lim.txt' # predictor_clean.py_bgl3_ss_v2_2880_lim
    # seq_and_score = 'bgl3_MLformat_not_ss_v2_2880_lim.txt' # predictor_clean.py_bgl3_not_ss_v2_2880_lim
    # seq_and_score = 'ube4b_MLformat_ss_v2_2880_lim.txt' # predictor_clean.py_ube4b_ss_v2_2880_lim
    # seq_and_score = 'ube4b_MLformat_not_ss_v2_2880_lim.txt' # predictor_clean.py_ube4b_not_ss_v2_2880_lim


    settings.seq = open(seq_and_score, 'r').readlines()[0].replace('\n', '').split(':')[0].split()

    main(seq_and_score, settings)