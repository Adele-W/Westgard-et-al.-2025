# Westgard-et-al.-2025
Contains scripts LabGrown and IDCrust. 

Author contact: adele.westgard@uit.no 

Link to article: TBD
The repository is created to host method developed for identifying laboratory-grown calcite (LabGrown) in _N.pachyderma_ and distinguish crust and lamellar calcite.

**LabGrown:**
Input: raw data from LA-ICP_MS as CSV files. 
Removes contaminated data, identifies and separates laboratory grown calcite based on the presence of a 135-Ba isotope spike. 
Output: csv files for each laser profile with clean, laboratory-grown data in ppm. 

**IDCrust:**
Input: 
  1) Output from LabGrown as CSV files. 
  2) Info spreadsheet (as xlsx) containing information about each laser ablation profile.
     Info in xlsx: laser shot ID, foraminifera ID, growth conditions (temperature, salinity, etc.) or core information (ID, depth, age, etc.)
Identifies and separates data in crust from lamellar calcite using breakpoints. Converts data from ppm to mmol/mol. 
Output: CSV file combining metadata (e.g., Foraminifera ID) and mean element/Ca in mmol/mol for crust, lamellar calcite and both components together.
        Can be modified to give full laser profile instead of mean. 

Citation: Westgård et al., 2025
          full citation TBC 
Copyright: Adele Westgård 2025. 
