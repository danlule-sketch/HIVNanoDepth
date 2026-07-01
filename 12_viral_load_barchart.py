#!/usr/bin/python

'''
This script reads viral load data from a CSV file and generates a bar chart of plasma 
viral load per sample. The y-axis is displayed on a logarithmic scale to accommodate 
wide dynamic ranges commonly observed in viral load measurements.
Adjust line 18 to reflect the path to the csv file
'''

# Import necessary libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

# Step 1: Load your CSV data into a DataFrame
data = pd.read_csv('/path_to_input/viral_load.csv')

# Step 2: Inspect the data
print(data.head())

# Step 3: Create the bar chart
plt.figure(figsize=(14, 6))  # Set figure size for better readability

# Create the bar chart using seaborn's barplot function
sns.barplot(x='IDNO', y='Plasma_Viral _Load(copies/mL)', data=data, color='steelblue')

# Step 4: Apply a logarithmic scale to the y-axis
plt.yscale('log')  # This will apply a log scale to the y-axis

# Step 5: Adjust x-axis labels to reduce crowding
tick_positions = np.arange(len(data))  # positions of bars
tick_labels = data['IDNO'].values    # original labels

# Show every 5th label
step = 5
plt.xticks(tick_positions[::step], tick_labels[::step], rotation=90)

# Step 6: Add titles and labels with bold text
plt.title('Viral Load by Sample', fontsize=20, fontweight='bold')
plt.xlabel('Sample', fontsize=16, fontweight='bold')
plt.ylabel('Viral Load (Log10 copies/mL)', fontsize=16, fontweight='bold')

# Step 7: Make axis values bold using FontProperties
font_props = font_manager.FontProperties(weight='bold')
plt.tick_params(axis='y', which='major', labelsize=15)
plt.yticks(fontproperties=font_props)

# Step 8: Save the plot to a file (e.g., PNG or PDF)
save_path = '/Users/daniellulebugembe/Documents/PhD/reading7_simulation/viral _loads/viral_load_barchart.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')  # Save with 300 dpi and tight bounding box

# Step 9: Show the plot
plt.tight_layout()  # Adjust layout to ensure labels fit properly
plt.show()
