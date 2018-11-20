# Contour : Pairwise Marine Distances

Minimum marine distances between sites in any region of the world.

## How it works?

A high-resolution polygon is converted to an infinite resistance surface. <br>
Minimum distances between sites are computed with a shortest path algorithm considering the infinite resistance of land and null resistance throughout the sea. <br>
The outcomes are a matrix of pairwise distances, a figure to visualize if sites are well represented in the study area and a figure depicting an example of a shortest marine distance. <br>
The main file with the sites should be structured as “Name Lon Lat” or “Name Lat Lon”. Coordinates must be in decimal degrees.

## Getting Started

Instructions to get the project up and running on your local machine.

### Prerequisites

Install the last verion of R available at [The Comprehensive R Archive Network](https://cran.r-project.org/).
Download a high resolution polygon depicting the surface of the world (e.g., Global Self-consistent Hierarchical High-resolution Shorelines; https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html) <br>

### Running the code

1. Open R and set the working directory (path to) <br>
2. Load the main function "contour" into memory (bellow) <br>
3. Run the  main function "contour"<br>
3.1 The main file with the locations should be text delimited<br>
3.2 Provide the path of the polygon depicting the surface of the world<br>
3.3 Define the delimiter and decimal character of the text file<br>
3.4 Define the main file structure: 1 to "Name Lon Lat" or 2 to "Name Lat Lon"<br>
3.5 Define if the text file has a header with the column names (TRUE or FALSE)<br>
3.6 Define the resolution of the study area and the buffer to use around the sites.  The buffer can be a simple value or a vector such as c(xmin,xmax,ymin,ymax).
3.7 Choose to export the results as a text delimited file (TRUE or FALSE) <br><br>

```
contour(  global.polygon = "Global_CostLine_HD_Polygon.shp" ,
file = "example.file.txt" , 
file.sep = "," ,
file.dec = "." ,
file.strucutre = 2 , 
file.header = FALSE ,
resolution = 0.01 ,
buffer = c(5,5,5,5) ,
export.file = TRUE   ) <br>
```

## Authors

* **Jorge Assis** @ [Bears Studio](https://www.bears.studio)

## License

Except where otherwise noted, the content on this repository is licensed under a [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).

## Citation

Assis, J., Castilho Coelho, N., Alberto, F., Valero, M., Raimondi, P., Reed, D., … Serrão, E. A. (2013). High and Distinct Range-Edge Genetic Diversity despite Local Bottlenecks. PLoS ONE, 8(7), e68646. https://doi.org/10.1371/journal.pone.0068646

