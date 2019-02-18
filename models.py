from keras.models import Model
from keras.layers import Input, Dense, Concatenate
import keras

def sorbnet(dim):
    X_input = Input(shape=(dim, ))
    X = Dense(16, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X_input) 
    X = Dense(8, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X)
    Y1 = Dense(4, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X)
    Y1 = Dense(4, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(Y1)
    Y1 = Dense(1, activation='sigmoid')(Y1)
    Y2 = Dense(4, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X)
    Y2 = Dense(4, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(Y2)
    Y2 = Dense(1, activation='sigmoid')(Y2)
    Y = Concatenate(1)([Y1, Y2])
    return Model(inputs=X_input, outputs=Y)

def shallow(dim):
    X_input = Input(shape=(dim, ))
    X = Dense(48, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X_input)
    Y1 = Dense(1, activation='sigmoid')(X)
    Y2 = Dense(1, activation='sigmoid')(X)
    Y = Concatenate(1)([Y1, Y2])
    return Model(inputs=X_input, outputs=Y)


def dense(dim):
    X_input = Input(shape=(dim, ))
    X = Dense(16, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X_input)
    X = Dense(8, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X)
    X = Dense(7, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X)
    X = Dense(6, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X) 
    Y1 = Dense(1, activation='sigmoid')(X)
    Y2 = Dense(1, activation='sigmoid')(X)
    Y = Concatenate(1)([Y1, Y2])
    return Model(inputs=X_input, outputs=Y)

def pvnet(dim):
    X_input = Input(shape=(dim, ))
    X = Dense(12, activation='elu', kernel_initializer='glorot_normal', bias_initializer='zeros')(X_input)
    X = Dense(1, activation='linear')(X)
    return Model(inputs=X_input, outputs=X)