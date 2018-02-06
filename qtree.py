#!/usr/bin/env python

import numpy as np

from osgeo import gdal, ogr, osr

class Qimage:
    """
    Quad image class containing the raster image within the quad. Also stores its corner points within the original
    image. The image is used to determine whether the quad is a leaf.
    """
    img = None  # numpy array
    transform = []

    def __init__(self, image, x1, x2, y1, y2):
        self.img = image
        self.x1 = x1; self.x2 = x2; self.y1 = y1; self.y2 = y2

    def set_val(self, val):
        self.val = val


class Quad:
    """
    Stores geographical vector information about the Qimage. Is used entirely for Shapefile outputs
    """
    ulx = 0; uly = 0; lrx = 0; lry = 0  # projected corner coordinates of the rectangle

    def __init__(self, ulx, uly, lrx, lry):
        self.ulx = ulx
        self.uly = uly
        self.lrx = lrx
        self.lry = lry

class Qtree:
    """
    Main class.
    """
    quad_list = []  # list of Quads making the resampled DEM

    def __init__(self, filename):
        """
        Initializes two images. One is the original, which stays the same throughout the program. The second is the
        quad-tree treated image.
        :param filename:
        """
        self.img = self.read_image(filename)
        self.img_q = self.img  # overwrites this for output .tif

    @staticmethod
    def read_image(filename):
        """
        Helper function to read a GeoTIFF image
        :param filename:
        :return:
        """
        ds = gdal.Open(filename)

        img = np.empty([ds.RasterXSize, ds.RasterYSize, ds.RasterCount])
        for i in range(ds.RasterCount):
            print "reading band " + str(i)
            img[:, :, i] = np.transpose(np.array(ds.GetRasterBand(i + 1).ReadAsArray()))

        img = Qimage(img[:, :, 0], 0, img.shape[0], 0, img.shape[1])  # currently only supports one band
        img.transform = ds.GetGeoTransform()

        ds = None

        return img


    def is_leaf(self, image):
        """
        Determines if a quad is a leaf or not based on some algebraic conditions. There are two conditions under which
        a quad is determined a leaf;
        (1) the size of the quad is less than a predefined value (in pixels)
        (2) the standard deviation of image values within the quad is less than a predefined value
        If either of these are true, then the function returns a true
        :param image:
        :return:
        """
        if image.img.shape[0] < self.n_minx or image.img.shape[1] < self.n_miny:
            image.set_val(np.mean(image.img))
            return True

        if np.std(image.img) < self.std_limit:
            image.set_val(np.mean(image.img))
            return True

        return False


    @staticmethod
    def divide(image):
        """
        Divides an image into four. These will be equal in size if the input image has side lengths divisible by 2.
        Otherwise they will be as close to equal sized as possible, meaning +1 pixel in side length
        :param image: an Qimage to be partioned
        :return: a list of partitioned images (Qimage)
        """
        x1 = image.img.shape[0] / 2
        x2 = image.img.shape[0]

        y1 = image.img.shape[1] / 2
        y2 = image.img.shape[1]

        # Return list of four squares starting quadrant 1, 2, 3, and then 4
        img1 = Qimage(image.img[ 0:x1,  0:y1], image.x1, image.x1 + x1, image.y1, image.y1 + y1)
        img2 = Qimage(image.img[ 0:x1, y1:y2], image.x1, image.x1 + x1, image.y1 + y1, image.y1 + y2)
        img3 = Qimage(image.img[x1:x2,  0:y1], image.x1 + x1, image.x1 + x2, image.y1, image.y1 + y1)
        img4 = Qimage(image.img[x1:x2, y1:y2], image.x1 + x1, image.x1 + x2, image.y1 + y1, image.y1 + y2)

        return [img1, img2, img3, img4]

    def quad_algorithm(self, img):
        for img in self.divide(img):
            if self.is_leaf(img):
                # add to image
                self.img_q.img[img.x1:img.x2, img.y1:img.y2] = img.val

                ulx = self.img.transform[0] + img.x1 * self.img.transform[1]
                uly = self.img.transform[3] + img.y1 * self.img.transform[5]
                lrx = self.img.transform[0] + img.x2 * self.img.transform[1]
                lry = self.img.transform[3] + img.y2 * self.img.transform[5]

                self.quad_list.append(Quad(ulx, uly, lrx, lry))
            else:
                self.quad_algorithm(img)

    def make(self, std_limit, n_minx, n_miny):
        """
        Make the quad-tree partitioning given the parameters
        :param std_limit: the standard deviation limit for a quad to be a leaf
        :param n_minx: minimum leaf side along x
        :param n_miny: minimum leaf sed along y
        :return:
        """
        self.n_minx = n_minx; self.n_miny = n_miny
        self.std_limit = std_limit
        self.quad_algorithm(self.img)

    def write_image(self, image):
        """
        Write the quad-tree partitioned image with mean value within quads
        :param image:
        :return:
        """
        src_ds = gdal.Open('./test/test.tif')

        drv = gdal.GetDriverByName('GTiff')
        dst_ds = drv.CreateCopy('./test/test_{}.tif'.format(self.std_limit), src_ds)

        dst_ds.GetRasterBand(1).WriteArray(np.transpose(image))  # write r-band to the    raster
        #ds.GetRasterBand(2).WriteArray(image)  # write g-band to the raster
        #ds.GetRasterBand(3).WriteArray(image)  # write b-band to the raster

        dst_ds.FlushCache()  # write to disk

        src_ds = None
        dst_ds = None

    def write_shp(self):
        """
        Write Shapefile to visualize quad-tree
        :return:
        """
        # set up the shapefile driver
        shp_driver = ogr.GetDriverByName("ESRI Shapefile")
        # create the data source

        shp_ds = shp_driver.CreateDataSource('./test/test{}_{}x{}_poly.shp'.format(self.std_limit, self.n_minx, self.n_miny))

        # create the spatial reference, WGS84
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(32632)

        # create and write the layer
        shp_lyr = shp_ds.CreateLayer('guide', srs, ogr.wkbPolygon)

        # Add an ID field
        field_id = ogr.FieldDefn("id", ogr.OFTInteger)
        shp_lyr.CreateField(field_id)

        id = 1
        for quad in self.quad_list:
            wkt = 'POLYGON (({} {}, {} {}, {} {}, {} {}, {} {}))'.format(
                quad.ulx, quad.uly,
                quad.lrx, quad.uly,
                quad.lrx, quad.lry,
                quad.ulx, quad.lry,
                quad.ulx, quad.uly  # Closes box (polygon)
            )

            geom = ogr.CreateGeometryFromWkt(wkt)

            feature = ogr.Feature(shp_lyr.GetLayerDefn())
            feature.SetGeometry(geom)
            feature.SetField("id", id)
            id += 1
            shp_lyr.CreateFeature(feature)
            feature = None

        shp_ds = None  # Save and close the data source

        # set up the shapefile driver
        shp_driver = ogr.GetDriverByName("ESRI Shapefile")
        # create the data source

        shp_ds = shp_driver.CreateDataSource('./test/test{}_{}x{}_point.shp'.format(self.std_limit, self.n_minx, self.n_miny))

        # create the spatial reference, WGS84
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(32632)

        # create and write the layer
        shp_lyr = shp_ds.CreateLayer('guide', srs, ogr.wkbPoint)

        # Add an ID field
        field_id = ogr.FieldDefn("id", ogr.OFTInteger)
        shp_lyr.CreateField(field_id)

        id = 1
        for quad in self.quad_list:
            wkt = 'POINT ({} {})'.format(
                (quad.ulx + quad.lrx) / 2,
                (quad.uly + quad.lry) / 2
            )

            geom = ogr.CreateGeometryFromWkt(wkt)

            feature = ogr.Feature(shp_lyr.GetLayerDefn())
            feature.SetGeometry(geom)
            feature.SetField("id", id)
            id += 1
            shp_lyr.CreateFeature(feature)
            feature = None

        shp_ds = None  # Save and close the data source



if __name__ == '__main__':
    tree = Qtree('./test/test.tif')
    tree.make(0.3, 16, 16)
    tree.write_image(tree.img_q.img)
    tree.write_shp()