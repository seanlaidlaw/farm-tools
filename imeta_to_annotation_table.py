#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
import re

# initialize dataframe we will append each sample's imeta to
imeta_dataframe = pd.DataFrame()

# for each argument given to script, convert to dictionary then add that to dataframe
for imeta_file in sys.argv[1:len(sys.argv)]:
    imeta_dict = {}

    # parse filename to remove directory name and imeta extension
    sample_id = re.sub(".*/", "", imeta_file)
    sample_id = re.sub("(\.imeta)?", "", sample_id)
    imeta_dict["seq_filename"] = sample_id


    # loop through imeta file assigning key,value pairs to each part
    with open(imeta_file) as f:
        imeta = f.readlines()

    attribute = ""
    value = ""
    for line in imeta:
        if line.startswith("attribute"):
            attribute = re.sub(r'^attribute: ', '', line.strip())
        elif line.startswith("value"):
            value = re.sub(r'^value: ', '', line.strip())
            imeta_dict[str(attribute)] = value

    # save key,value pairs to pandas dataframe
    imeta_dataframe = imeta_dataframe.append(imeta_dict, ignore_index=True)

# set rewnames of dataframe to be the bam filename
imeta_dataframe = imeta_dataframe.set_index("seq_filename")

# save dataframe as tsv
imeta_dataframe.to_csv("imeta_information.tsv", sep='\t')

imeta_dataframe=imeta_dataframe[[i for i in imeta_dataframe if len(set(imeta_dataframe[i]))>1]]

# save dataframe as tsv
imeta_dataframe.to_csv("imeta_information_nosame.tsv", sep='\t')
