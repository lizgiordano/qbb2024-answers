#!/usr/bin/env python


import numpy as np
import scipy
import matplotlib.pyplot as plt # Plotting package for python 
import imageio
import imageio.v3
import plotly.express as px 
import plotly

# transform array dimension (axis) - to reshape array to make it more manageable
na = np.newaxis # Numpy constant

# Exercise 1 Load images, use a for loop for going through each dimension
# load images into memory, process them, and store them in a dictionary.
# Read image files based on gene names, field , and channels, normalizes the image data, and stores as NumPy arrays.

# creating a dictionary to make it easier to group the 3 dimensions: gene name, field, channel color
images = {}
image_dir = "/Users/cmdb/qbb2024-answers/qb24/week10/"
genes = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
channels = ["DAPI", "PCNA", "nascentRNA"]
# for loop to go through each genename, field, and channel
for gene in genes:  # loops through each gene
    for field in fields:  # for each gene, loops through field 
        image = np.zeros((520, 616, 3), np.uint16)  # create array with those dimensions
        for i, channel in enumerate(channels):
            image_name = image_dir + gene + '_' + field + '_' + channel + '.tif'
            try:
                image[:, :, i] = imageio.v3.imread(image_name)
            except FileNotFoundError:
                print(f"Warning: File {image_name} not found.")
                continue
        image = image.astype(np.float32)  # resize
        for i in range(3):
            image[:, :, i] -= np.amin(image[:, :, i])
            image[:, :, i] /= np.amax(image[:, :, i])
        # Rescale and clip values between 0 and 255
        image = np.clip(image * 255, 0, 255).astype(np.uint8)
        images[f"{gene}_{field}"] = image
    

# exercise 2
# Step 2.1 For each image, create a binary mask from the DAPI channel

# mask differentiates between nucleus and cytoplasm
# gotta find nucleus by DAPI
dapi_mean = np.mean(image[:, :, 0]) 
mask = image[:, :, 0] >= dapi_mean



# step 2.2
# Find labels for each image based on the DAPI mask from step 2.1
# should get an array of shape X.Y
#zeros in all non-cell positions and positive numbers denoting pixels belonging to a given nucleus
# go over the binary mask, pixel by pixel, to assign labels to connected regions.

# manually assign labels pixel by pixel so they can connect to the same region
# goes through neighboring pixels (up, left, right) to make sure connected pixels belong to the same region

def find_labels(mask): # assign labels, make a list
    l = 0
    labels = np.zeros(mask.shape, np.int32)
    equivalence = [0]
    if mask[0, 0]: # gotta go through each region, up, right, left
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    for x in range(1, mask.shape[0]): # for nonzero rows
        if mask[x, 0]:
            if mask[x - 1, 0]:
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    if mask[x - 1, y - 1]:
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        if mask[x, -1]:  # Check the last pixel 
            if mask[x - 1, -1]:
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                labels[x, -1] = equivalence[labels[x, -2]]
            else:  # add new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = np.array(equivalence)
    for i in range(1, len(equivalence))[::-1]:
        labels[np.where(labels == i)] = equivalence[i]
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        labels[np.where(labels == j)] = i
    return labels

labels = find_labels(mask)

# plt.imshow(labels)
#plt.show()



# Step 2.3
# Filter out labeled outliers based on size
# use the filter_by_size function from the live coding for this step.

# set minsize and maxsize to filter by pixel size
size = np.bincount(labels.ravel())
for i in range(1, np.amax(labels)+1):
    where = np.where(labels == i)
    if size[i] < 100:
        labels[where] = 0
# use filter by size and live coding   
#  leave outliers as background     
def filter_by_size(labels, minsize, maxsize):
    sizes = np.bincount(labels.ravel())
    for i in range(1, sizes.shape[0]):
        if sizes[i] < minsize or sizes[i] > maxsize:
            where = np.where(labels == i)
            labels[where] = 0
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        labels[np.where(labels == j)] = i
    return labels


# Exercise 3: Score the PCNA and nascent RNA signal in each nucleus and plot them

# step 3.1 
# Find the mean signal for each nucleus from the PCNA and nascent RNA channels

# go through genes and fields, and use DAPI mask to get nuclei
# then use the filter by size so pixels are grouped into nuclei
# then loop through each labelled nuclei
# change to PCNA 
# get ratio of nascent to pcna
# save file with labelled tabs
mean = [] 
for gene in genes: 
    for field in fields: 
        image = images[f"{gene}_{field}"] 
        dapi_mask = mask  
        labels = find_labels(dapi_mask)  
        labels = filter_by_size(labels, minsize=100, maxsize=10000)  
        mean_nascentRNA = []
        mean_pcna = []
        ratios = []
        for i in range(1, np.amax(labels) + 1): 
            where = np.where(labels == i) 
            mean_nascentRNA.append(np.mean(image[where[0], where[1], 1])) 
            mean_pcna.append(np.mean(image[where[0], where[1], 2])) 
        for nRNA, pcna in zip(mean_nascentRNA, mean_pcna): 
            if pcna != 0:  
                ratios.append(np.log2(nRNA / pcna)) 
            else:
                ratios.append(0)  
        for mnRNA, mPCNA, ratio in zip(mean_nascentRNA, mean_pcna, ratios): 
            mean.append([gene, mnRNA, mPCNA, ratio]) 
output_filename = "mean_mucleus_signal.txt" 
with open(output_filename, 'w') as file:
    file.write("Gene\tNascentRNA\PCNA\tRatio\n") 
    for mean in means: 
        gene = mean[0]  
        mean_nascentRNA = mean[1] 
        mean_pcna = mean[2] 
        ratio = mean[3] 
        file.write(f"{gene_field}\t{mean_nascentRNA:.4f}\t{mean_pcna:.4f}\t{ratio:.4f}\n")
        


