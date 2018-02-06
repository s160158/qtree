# qtree
A tool to generate a quad-tree on a GeoTIFF spatial image. The idea is to obtain points representing the change in the image to constitue vertices in a Delayney triangulation or assist the MIKE Mesh Generator.

A total of three products are produced; quad center points, quad outlines, and a GeoTIFF containing the mean values within the quads. 

The quad center point shapefile can be converted with my shp2mdf to produce a file compatible with MIKE Mesh Generator. An outline of the area will be needed to make the mesh.




