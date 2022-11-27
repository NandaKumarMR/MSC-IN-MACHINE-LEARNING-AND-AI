import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import joblib

# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/SubstructureFingerprintCount.xml -dir ./ -file descriptors_output_r.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output_c.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data_r,input_data_c):
    # Reads in saved regression model
    load_model_r = pickle.load(open('acetylcholinesterase_model_regressor.pkl', 'rb'))
    load_model_c = pickle.load(open('acetylcholinesterase_model_classifier.pkl', 'rb'))
    # Apply model to make predictions
    prediction_r = load_model_r.predict(input_data_r)
    prediction_c = load_model_c.predict(input_data_c)
    st.header('**Prediction output**')
    prediction_output_r = pd.Series(prediction_r, name='pIC50')
    for i in df.pIC50:
        if float(i) <= 5.0:
            bioactivity_threshold.append("Inactive")
        elif float(i) >= 6.0:
            bioactivity_threshold.append("Active")
        else:
            bioactivity_threshold.append("Intermediate")
    bioactivity_class_r = pd.Series(bioactivity_threshold, name = 'Bioactivity Class Regressor')
    prediction_output_c = pd.Series(prediction_c, name=''Bioactivity Class Classifier'')
    molecule_name = pd.Series(load_data[1], name='Molecule Name')
    df = pd.concat([molecule_name, prediction_output_r,bioactivity_class_r, prediction_output_c], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# Logo image
image = Image.open('App_logo.png')

st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Bioactivity Prediction App (Acetylcholinesterase)

This app allows you to predict the bioactivity towards inhibting the `Acetylcholinesterase` enzyme. `Acetylcholinesterase` is a drug target for Alzheimer's disease.

**Credits**
- App built in `Python` + `Streamlit`
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
---
""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/NandaKumarMR/MSC-IN-MACHINE-LEARNING-AND-AI/main/Input_molecules.txt)
""")

if st.sidebar.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors Regression**')
    desc_r = pd.read_csv('descriptors_output_r.csv')
    st.write(desc_r)
    st.write(desc_r.shape)
    
    st.header('**Calculated molecular descriptors Classification**')
    desc_c = pd.read_csv('descriptors_output_c.csv')
    st.write(desc_c)
    st.write(desc_c.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built regression models**')
    Xlist = list(pd.read_csv('descriptor_list_r.csv').columns)
    desc_subset_r = desc_r[Xlist]
    st.write(desc_subset_r)
    st.write(desc_subset_r.shape)
    
    st.header('**Subset of descriptors from previously built classification models**')
    Xlist = list(pd.read_csv('descriptor_list_c.csv').columns)
    desc_subset_c = desc_c[Xlist]
    st.write(desc_subset_c)
    st.write(desc_subset_c.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset_r,desc_subset_c)
else:
    st.info('Upload input data in the sidebar to start!')
