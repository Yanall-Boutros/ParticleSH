# ---------------------------------------------------------------------
# Import Statements
# ---------------------------------------------------------------------
import tensorflow as tf
import numpy as np
import pickle as pic
import os
import os.path
import time
# ---------------------------------------------------------------------
# Function Definitions
# ---------------------------------------------------------------------
def pred_comp(pred, real):
    # Predictions compare takes a predictions array as input,
    # and compares its results to an expected output (real array)
    # Returns the percentage of predictions which were accurate.
    c = 0
    predictions = np.round(pred)
    elems = len(predictions)
    for i in range(elems):
        if predictions[i] == real[i]:
            c += 1
    return c / elems

def get_care_packages():
    base_name = "care_package"
    i = 0
    # If there are no care packages, then wait
    while not os.path.isfile(base_name+str(i)):
        print("No care_package file detected, waiting 10 seconds...")
        time.sleep(10)
        # quick check to see if it's time to stop
        if not np.load(open("control"), "rb"): exit()
    # iterate through all the care_packagen files where n is an arbitrary
    # number for valid carepackages, then export all as a list, deleting as
    # we append to the list
    care_pallet = [] # A pallet contains many packages
    while os.path.isfile(base_name+str(i)):
        care_pallet.append(np.load(open(base_name+str(i), "rb")))
        os.remove(base_name+str(i))
        i += 1
    return care_pallet
# ---------------------------------------------------------------------
# Check if model exists, if so, load it. Otherwise, make new network
# ---------------------------------------------------------------------
if os.path.isfile("ff_model"):
    ffmodel = tf.keras.models.load_model("ffmodel")
else:
    ffmodel = tf.keras.models.Sequential()
    ffmodel.add(tf.keras.layers.Flatten())
    ffmodel.add(tf.keras.layers.Dense(64, activation=tf.nn.relu))
    ffmodel.add(tf.keras.layers.Dense(25, activation=tf.nn.relu))
    ffmodel.add(tf.keras.layers.Dense(1, activation=tf.nn.sigmoid))

ffmodel.compile(loss='binary_crossentropy',
                optimizer='Adam',
                metric=['acurracy'])
# -----------------------------------------------------------------------
# Train and Test the Network
# -----------------------------------------------------------------------
while np.load(open("control", "rb")):
    care_pallet = get_care_packages()
    for care_package in care_pallet:
        T_i, T_o = care_package[0], care_package[1] 
        Test_i, Test_o = care_package[2], carepackage[3]
        model.fit(T_i, T_o, epochs=100)
        model.save("first_model")
        predicitons = model.predict(Test_i)
        print("Neural Network correctly evaluated ", 100*pred_comp(predicitons, Test_o)
              ,"% of test data")
print("Training Halted")
