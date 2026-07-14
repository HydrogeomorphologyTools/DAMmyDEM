# DAMmyDEM

**DAMmyDEM** is a hydrogeomorphological Python script designed to generate a "virtual dam" and its corresponding water reservoir directly within a Digital Elevation Model (DEM). 

## 📖 How it Works
The script modifies a baseline DEM to simulate the topography of a basin retained by an artificial dam. It requires an input DEM of the study area and a polygonal shapefile describing the dam's footprint (e.g., a valley-transverse rectangular polygon).

### Core Workflow:
1. **Elevation Extraction:** The script masks the input DEM using the dam shapefile to extract the terrain pixels falling inside the polygon.
2. **Crest Level Calculation:** It identifies the highest elevation value among these extracted pixels. This maximum elevation acts as the dam crest and defines the maximum water level of the reservoir.
3. **Water Backpropagation (TauDEM):** Taking the calculated maximum water level as a reference, the script utilizes **TauDEM** algorithms to backpropagate a perfectly flat water surface upstream along the drainage network. 
4. **DEM Modification:** The script outputs an updated DEM where the original valley topography upstream of the dam is mathematically "filled" up to the newly calculated water level.

## 🛠️ System Dependencies
The script's core hydrological routing and filling operations require **TauDEM** (Terrain Analysis Using Digital Elevation Models). 

* **TauDEM**: Must be installed and accessible on your system. It is available at the [official website](https://hydrology.usu.edu/taudem/taudem5/) or via [GitHub Releases](https://github.com/dtarb/TauDEM/releases).
* **Python**: The tool is written in Python. Ensure your environment has the necessary geospatial libraries installed (such as `rasterio` and `geopandas`) to handle the raster and vector inputs.

## 🗂️ Data Preparation & Usage
To ensure accurate calculations and prevent artifacts in the backpropagation, your input data must be properly formatted:

* **DEM Filtering:** The input `.tif` files must be pre-processed so that only valid elevation values strictly greater than 0 and lower than 5000 are considered. Any other values outside this range must be explicitly converted to `NoData` before running the tool.
* **Required Inputs:**
  * The filtered **DEM** (Raster `.tif` format).
  * The **Dam Shapefile** (Polygonal vector format).
* **Execution:** All file paths, file names, and output directories must be specified inside the `param_dam_dem.json` configuration file before launching the script.

## 👥 Credits
Developed by **Stefano Crema** and **Marco Cavalli** (2024).

## 📄 License
This project is licensed under the **GNU General Public License v2.0 (GPL-2.0)**.
