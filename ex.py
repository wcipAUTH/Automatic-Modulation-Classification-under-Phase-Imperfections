import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import h5py
import mat73
from tensorflow import keras
from keras.models import Sequential
from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense, Dropout
from keras.utils import np_utils
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score
import time
import sys

# Save original sys.stdout
original_stdout = sys.stdout

# Open file to save output

    # Redirect sys.stdout to the file
   
beaft = 'after'
k = 10
print(k)
print(beaft)
nn_type = 'LeNet'
print(nn_type)
# Load the training, validation, and test data from .mat files
train_data = mat73.loadmat('datasets/datasets_train_hqam_k={}.mat'.format(k))
X_train = train_data['Images_total_{}'.format(beaft)]
y_train = train_data['labels_onehot_total_{}'.format(beaft)]

val_data = mat73.loadmat('datasets/datasets_val_hqam_k={}.mat'.format(k))
X_val = val_data['Images_total_{}'.format(beaft)]
y_val = val_data['labels_onehot_total_{}'.format(beaft)]

test_data = mat73.loadmat('datasets/datasets_test_hqam_k={}.mat'.format(k))
X_test = test_data['Images_total_{}'.format(beaft)]
y_test = test_data['labels_onehot_total_{}'.format(beaft)]
z_test = test_data['Image_order_total_{}'.format(beaft)]

    #change the number of images accordingly
X_train = np.reshape(X_train, (448*16, 256, 256, 1))
X_val = np.reshape(X_val, (224*16 ,  256, 256, 1))
X_test = np.reshape(X_test, (1120*16,256, 256, 1))


    # Define the AlexNet-5 model architecture
if (nn_type=='AlexNet-asymmetric'):
    model = Sequential([
            Conv2D(filters=6, kernel_size=(5, 5), activation='relu', input_shape=(256, 256, 1)),
            MaxPooling2D(pool_size=(2, 2)),
            Conv2D(filters=16, kernel_size=(3, 6), activation='relu'),
            MaxPooling2D(pool_size=(2, 2)),
            Conv2D(filters=16, kernel_size=(3, 9), activation='relu'),
            Conv2D(filters=16, kernel_size=(3, 3), activation='relu'),
            Conv2D(filters=16, kernel_size=(3, 3), activation='relu'),
            MaxPooling2D(pool_size=(2, 2)),
            Flatten(),
            Dense(units=120, activation='relu'),
            Dense(units=84, activation='relu'),
            Dense(units=16,activation='softmax')

         ])
elif (nn_type == 'AlexNet'):
    model = Sequential([
            Conv2D(filters=6, kernel_size=(5, 5), activation='relu', input_shape=(256, 256, 1)),
            MaxPooling2D(pool_size=(2, 2)),
            Conv2D(filters=16, kernel_size=(3, 3), activation='relu'),
            MaxPooling2D(pool_size=(2, 2)),
            Conv2D(filters=16, kernel_size=(3, 3), activation='relu'),
            Conv2D(filters=16, kernel_size=(3, 3), activation='relu'),
            Conv2D(filters=16, kernel_size=(3, 3), activation='relu'),
            MaxPooling2D(pool_size=(2, 2)),
            Flatten(),
            Dense(units=120, activation='relu'),

            Dense(units=84, activation='relu'),
            Dense(units=16,activation='softmax')

        ])
elif (nn_type == 'LeNet'):
    model = Sequential([
            Conv2D(filters=6, kernel_size=(5, 5), activation='relu', input_shape=(256, 256, 1)),
            MaxPooling2D(pool_size=(2, 2)),
            Conv2D(filters=16, kernel_size=(3, 3), activation='relu'),
            MaxPooling2D(pool_size=(2, 2)),
            Flatten(),
            Dense(units=120, activation='relu'),
            Dense(units=84, activation='relu'),
            Dense(units=16, activation='softmax')
        ])
elif (nn_type== 'LeNet-asymmetric'):
    model = Sequential([
            Conv2D(filters=6, kernel_size=(5, 5), activation='relu', input_shape=(256, 256, 1)),
            MaxPooling2D(pool_size=(2, 2)),
            Conv2D(filters=16, kernel_size=(3, 6), activation='relu'),
            MaxPooling2D(pool_size=(2, 2)),
            Flatten(),
            Dense(units=120, activation='relu'),
            Dense(units=84, activation='relu'),
            Dense(units=16, activation='softmax')
        ])
else:
    print('wrong input')
    

'''
model = Sequential([
        Conv2D(filters=96, kernel_size=(5, 5), strides = (2,2), padding="valid", activation='relu', input_shape=(256, 256, 1)),
        MaxPooling2D(pool_size=(3, 3),strides = (2,2)),
        Conv2D(filters=256, kernel_size=(5, 5), padding="same",activation='relu'),
        MaxPooling2D(pool_size=(3, 3), strides= (2,2)),
        Conv2D(filters=384, kernel_size=(3, 3), padding="same", activation='relu'),
        Conv2D(filters=384, kernel_size=(3, 3), padding="same", activation='relu'),
        Conv2D(filters=256, kernel_size=(3, 3), padding="same", activation='relu'),
        MaxPooling2D(pool_size=(4, 4),strides= (2,2)),
        Flatten(),
        Dense(units=512, activation='relu'),
        Dropout(0.5),
        Dense(units=512, activation='relu'),
        Dropout(0.5),
        Dense(units=16, activation='relu')
    ])
'''

    # Compile the LeNet model
model.compile(keras.optimizers.Adam(learning_rate=1e-4),loss='categorical_crossentropy', metrics=['accuracy'])

    # summarize the model
model.summary()

    # Train the model
history = model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=1, batch_size=64)
'''
test_loss, test_acc = model.evaluate(X_test,y_test)
'''
av_time = 0
reps = 3
for i in range(1,reps):
    start = time.time()
    test_loss, test_acc = model.evaluate(X_test,y_test)
    end = time.time()
    av_time_temp = end-start
    av_time = av_time + av_time_temp

print(av_time/reps)


from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import numpy as np

    # Get predicted labels for test set
y_pred = model.predict(X_test)

    # Convert predicted and true labels to one-dimensional arrays
y_pred_labels = np.argmax(y_pred, axis=1)
y_true_labels = np.argmax(y_test, axis=1)

    # Create confusion matrix
cm = confusion_matrix(y_true_labels, y_pred_labels)

    # Compute the classification report
cr = classification_report(y_true_labels, y_pred_labels)

    # Calculate success percentage for each case in the confusion matrix
cm_percentage = np.round(cm.astype('float') / cm.sum(axis=1)[:, np.newaxis], decimals=2)

    # Compute the accuracy score
acc = accuracy_score(y_true_labels, y_pred_labels)

print("Confusion Matrix:\n", cm)
print("Classification Report:\n", cr)
print("Overall Accuracy Score:", acc)

    # Print confusion matrix as an image
fig, ax = plt.subplots(figsize=(10, 7))
im = ax.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
ax.figure.colorbar(im, ax=ax)

    # Set axis labels
    # change the name of the labels accordingly
classes = ['4-ASK','8-ASK','BPSK','QPSK','4-HQAM','16-HQAM','64-HQAM','16-QAM','32-QAM','64-QAM','128-QAM','256-QAM','16-APSK','32-APSK','64-APSK','128-APSK']
ax.set(xticks=np.arange(cm.shape[1]),
        yticks=np.arange(cm.shape[0]),
        xticklabels=classes, yticklabels=classes,
        xlabel='Predicted label',
        ylabel='True label')

    # Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        ax.text(j, i, "{:.0%}".format(cm_percentage[i, j]),  # Display success percentage
                    ha="center", va="center",
                    color="white" if cm[i, j] > cm.max() / 2. else "black")
        fig.tight_layout()

# Save the confusion matrix plot as an EPS file
plt.savefig('confusion_matrix_after_net1.eps', format='eps')

    # Automatically adjust subplot parameters to fit the plot within the figure area
plt.tight_layout()

    # Display the plot
plt.show()


    #change snr range accordingly
    #SNR_dB = [-25,-20,-15,-10,0,5]

SNR_dB = [0,5,10,15,20,25,30]
acc = []
acc_X_test = []
acc_Y_test = []

for i in range(len(SNR_dB)):
    acc_snr = SNR_dB[i]
    idx_acc_snr = []

    for k in range(len(z_test)):
        if int(z_test[k]) == int(acc_snr):
            idx_acc_snr.append(k)
    acc_X_test = X_test[idx_acc_snr,:,:,:]
    acc_Y_test = y_test[idx_acc_snr,:]

            # evaluate the model
    _, accuracy_snr = model.evaluate(acc_X_test, acc_Y_test, batch_size=66)
    print('SNR ' + str(acc_snr) + 'dB: ' + str(accuracy_snr))

        # Add the accuracy to the list
    acc.append(accuracy_snr)

    # Plot training and validation accuracy values
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('Model accuracy')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.legend(['Train', 'Validation'], loc='upper left')
plt.show()

plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Validation'], loc='upper left')
plt.show()

plt.plot(SNR_dB, acc)
#change snr range accordingly
plt.xlim([0,30])
plt.ylim([min(acc)-0.01 ,  max(acc)+0.01])
plt.xlabel('SNR (dB)')
plt.ylabel('Test accuracy')
plt.title('Test accuracy vs. SNR')
plt.show()

# Create a list of (x, y) coordinates
data_points = list(zip(SNR_dB, acc))

# Define the output file path
output_file = 'test_accuracy_vs_SNR_after_net1.dat'

# Save the data points to the .dat file
with open(output_file, 'w') as f:
    for point in data_points:
        f.write(f'{point[0]}\t{point[1]}\n')

# Restore sys.stdout to original value
sys.stdout = original_stdout
