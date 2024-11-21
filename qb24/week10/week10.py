!/usr/bin/env python


import numpy
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly
import plotly.express as px

# transform array dimension (axis) - to reshape array to make it more manageable
na = np.newaxis # Numpy constant


# Exercise 1 Load images, use a for loop for going through each dimension

# creating a dictionary to make it easier to group the 3 dimensions: gene name, field, channel color
images = {}

# now name the dimensions/variables from the file name that I want for the for loop 
names = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
channels = ["DAPI", "PCNA", "nascentRNA"]

# make a for loop that goes through the name, field, and channel
for name in names: # take the gene name from the file
    images[name] = {} # store it
    for field in fields:
        images[name][field] = {}  
        for channel in channels:
            # put them together for the filename
            filepath = f"~/qbb2024-answers/qb24/week10/{name}_{field}_{channel}"
            try:
                img = imageio.v3.imread(filepath).astype(np.uint16) # Read the image and store it in the dictionary
                images[f"{name}_{field}_{channel}"] = img
            except FileNotFoundError:
                # exclude missing files 
                images[gene][field][channel] = None
                #print(f"File not found: {filename}")

# need another for loop to make a list of 3D image arrays
all_img_arrays = [] # makes list

# Loop through genes, fields, channels
for name in names:
    for field in fields:
        # make array to store channels
        img_array = np.zeros((test_img.shape[0], test_img.shape[1], 3), np.uint16)
        for i, channel in enumerate(channels):
            try:
                # Retrieve image
                img = images[f"{name}_{field}_{channel}"]
                if img is None:
                    raise KeyError(f"Image missing: {name}_{field}_{channel}")
                # process and store in array
                temp = img.astype(np.float32)
                temp -= np.amin(temp)
                temp /= np.amax(temp)
                img_array[:, :, i] = (temp * 65535).astype(np.uint16)  
            except KeyError as e:
                # print(e)
        # add the new processed array to the list
        all_img_arrays.append(img_array)

# Convert to NumPy array
all_img_arrays = np.array(all_img_arrays)

# exercise 2
# Step 2.1 For each image, create a binary mask from the DAPI channel

# mask differentiates between nucleus and cytoplasm
# gotta find nucleus by DAPI
mask = []
for img in all_img_arrays:
    DAPI_channel = img[:, :, 0]
    DAPI_mean = np.mean(DAPI_channel)
    DAPI_mask = DAPI_channel >= DAPI_mean
    mask.append(DAPI_mask)
plt.imshow(mask[0])
plt.show()








