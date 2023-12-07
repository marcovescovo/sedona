/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

package org.apache.sedona.core.showcase;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.sedona.common.enums.FileDataSplitter;
import org.apache.sedona.core.enums.GridType;
import org.apache.sedona.core.enums.IndexType;
import org.apache.sedona.core.formatMapper.shapefileParser.ShapefileRDD;
import org.apache.sedona.core.serde.SedonaKryoRegistrator;
import org.apache.sedona.core.spatialOperator.JoinQuery;
import org.apache.sedona.core.spatialOperator.KNNQuery;
import org.apache.sedona.core.spatialOperator.RangeQuery;
import org.apache.sedona.core.spatialRDD.CircleRDD;
import org.apache.sedona.core.spatialRDD.PointRDD;
import org.apache.sedona.core.spatialRDD.PolygonRDD;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.storage.StorageLevel;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

// TODO: Auto-generated Javadoc

/**
 * The Class Example.
 */
public class Example
        implements Serializable
{

    /**
     * The sc.
     */
    public static JavaSparkContext sc;

    /**
     * The geometry factory.
     */
    static GeometryFactory geometryFactory;

    /**
     * The Point RDD input location.
     */
    static String PointRDDInputLocation;

    /**
     * The Point RDD offset.
     */
    static Integer PointRDDOffset;

    /**
     * The Point RDD num partitions.
     */
    static Integer PointRDDNumPartitions;

    /**
     * The Point RDD splitter.
     */
    static FileDataSplitter PointRDDSplitter;

    /**
     * The Point RDD index type.
     */
    static IndexType PointRDDIndexType;

    /**
     * The object RDD.
     */
    static PointRDD objectRDD;

    /**
     * The Polygon RDD input location.
     */
    static String PolygonRDDInputLocation;

    /**
     * The Polygon RDD start offset.
     */
    static Integer PolygonRDDStartOffset;

    /**
     * The Polygon RDD end offset.
     */
    static Integer PolygonRDDEndOffset;

    /**
     * The Polygon RDD num partitions.
     */
    static Integer PolygonRDDNumPartitions;

    /**
     * The Polygon RDD splitter.
     */
    static FileDataSplitter PolygonRDDSplitter;

    /**
     * The query window RDD.
     */
    static PolygonRDD queryWindowRDD;

    /**
     * The join query partitioning type.
     */
    static GridType joinQueryPartitioningType;

    /**
     * The each query loop times.
     */
    static int eachQueryLoopTimes;

    /**
     * The k NN query point.
     */
    static Point kNNQueryPoint;

    /**
     * The range query window.
     */
    static Envelope rangeQueryWindow;

    static String ShapeFileInputLocation;

    /**
     * The main method.
     *
     * @param args the arguments
     */
    public static void main(String[] args)
    {
        SparkConf conf = new SparkConf().setAppName("SedonaRunnableExample").setMaster("local[2]");
        conf.set("spark.serializer", KryoSerializer.class.getName());
        conf.set("spark.kryo.registrator", SedonaKryoRegistrator.class.getName());

        sc = new JavaSparkContext(conf);
        Logger.getLogger("org").setLevel(Level.WARN);
        Logger.getLogger("akka").setLevel(Level.WARN);

        String resourceFolder = System.getProperty("user.dir") + "/src/test/resources/";

        PointRDDInputLocation = resourceFolder + "arealm-small.csv";
        PointRDDSplitter = FileDataSplitter.CSV;
        PointRDDIndexType = IndexType.RTREE;
        PointRDDNumPartitions = 5;
        PointRDDOffset = 0;

        PolygonRDDInputLocation = resourceFolder + "primaryroads-polygon.csv";
        PolygonRDDSplitter = FileDataSplitter.CSV;
        PolygonRDDNumPartitions = 5;
        PolygonRDDStartOffset = 0;
        PolygonRDDEndOffset = 8;

        geometryFactory = new GeometryFactory();
        kNNQueryPoint = geometryFactory.createPoint(new Coordinate(500000, 500000));
        rangeQueryWindow = new Envelope(-90.01, -80.01, 30.01, 40.01);
        joinQueryPartitioningType = GridType.QUADTREE;
        eachQueryLoopTimes = 5;

        ShapeFileInputLocation = resourceFolder + "shapefiles/polygon";

        try {
            /*testSpatialRangeQuery();
            testSpatialRangeQueryUsingIndex();*/
            testSpatialKnnQuery();
            /*testSpatialKnnQueryUsingIndex();
            testSpatialJoinQuery();
            testSpatialJoinQueryUsingIndex();
            testDistanceJoinQuery();
            testDistanceJoinQueryUsingIndex();
            testLoadShapefileIntoPolygonRDD();*/
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println("DEMOs failed!");
            return;
        }
        sc.stop();
        System.out.println("All DEMOs passed!");
    }

    /**
     * Test spatial range query.
     *
     * @throws Exception the exception
     */
    public static void testSpatialRangeQuery()
            throws Exception
    {
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);
        objectRDD.rawSpatialRDD.persist(StorageLevel.MEMORY_ONLY());
        for (int i = 0; i < eachQueryLoopTimes; i++) {
            long resultSize = RangeQuery.SpatialRangeQuery(objectRDD, rangeQueryWindow, false, false).count();
            assert resultSize > -1;
        }
    }

    /**
     * Test spatial range query using index.
     *
     * @throws Exception the exception
     */
    public static void testSpatialRangeQueryUsingIndex()
            throws Exception
    {
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);
        objectRDD.buildIndex(PointRDDIndexType, false);
        objectRDD.indexedRawRDD.persist(StorageLevel.MEMORY_ONLY());
        for (int i = 0; i < eachQueryLoopTimes; i++) {
            long resultSize = RangeQuery.SpatialRangeQuery(objectRDD, rangeQueryWindow, false, true).count();
            assert resultSize > -1;
        }
    }

    /**
     * Test spatial knn query.
     *
     * @throws Exception the exception
     */
    @SuppressWarnings("null")
	public static void testSpatialKnnQuery()
            throws Exception
    {
    	////////////////////////////////////////////////////////////////
        // Define a formatter for the timestamp (optional)
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss.SSS");
        ////////////////////////////////////////////////////////////////
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);
        objectRDD.rawSpatialRDD.persist(StorageLevel.MEMORY_ONLY());
        ArrayList<Point> result_final = new ArrayList<Point>();
        //for (int i = 0; i < eachQueryLoopTimes; i++) {
        //for (int i = 0; i < objectRDD.approximateTotalCount; i++) {
        	//List<Point> result = KNNQuery.SpatialKnnQuery(objectRDD, kNNQueryPoint, 1000, false);
        LocalDateTime startDateTime = LocalDateTime.now();
        String formattedTimestamp = startDateTime.format(formatter);
        System.out.println("Start Timestamp: " + formattedTimestamp);
    	List<Point> result = KNNQuery.SpatialKnnQuery(objectRDD, kNNQueryPoint, 100000, false);
    	LocalDateTime endDateTime = LocalDateTime.now();
        formattedTimestamp = endDateTime.format(formatter);
        System.out.println("End Timestamp: " + formattedTimestamp);
        assert result.size() > -1;
        if (result.size() > 0) {
        	for (int j = 0; j < result.size(); j++) {
        		if (result_final.contains(result.get(j)) == false) {}
                result_final.add(result.get(j));	
        	}	
        }
        
        System.out.println(result_final.size());
        System.out.println();
    }

    /**
     * Test spatial knn query using index.
     *
     * @throws Exception the exception
     */
    public static void testSpatialKnnQueryUsingIndex()
            throws Exception
    {
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);
        objectRDD.buildIndex(PointRDDIndexType, false);
        objectRDD.indexedRawRDD.persist(StorageLevel.MEMORY_ONLY());
        for (int i = 0; i < eachQueryLoopTimes; i++) {
            List<Point> result = KNNQuery.SpatialKnnQuery(objectRDD, kNNQueryPoint, 1000, true);
            assert result.size() > -1;
        }
    }

    /**
     * Test spatial join query.
     *
     * @throws Exception the exception
     */
    public static void testSpatialJoinQuery()
            throws Exception
    {
        queryWindowRDD = new PolygonRDD(sc, PolygonRDDInputLocation, PolygonRDDStartOffset, PolygonRDDEndOffset, PolygonRDDSplitter, true);
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);

        objectRDD.spatialPartitioning(joinQueryPartitioningType);
        queryWindowRDD.spatialPartitioning(objectRDD.getPartitioner());

        objectRDD.spatialPartitionedRDD.persist(StorageLevel.MEMORY_ONLY());
        queryWindowRDD.spatialPartitionedRDD.persist(StorageLevel.MEMORY_ONLY());
        for (int i = 0; i < eachQueryLoopTimes; i++) {
            long resultSize = JoinQuery.SpatialJoinQuery(objectRDD, queryWindowRDD, false, true).count();
            assert resultSize > 0;
        }
    }

    /**
     * Test spatial join query using index.
     *
     * @throws Exception the exception
     */
    public static void testSpatialJoinQueryUsingIndex()
            throws Exception
    {
        queryWindowRDD = new PolygonRDD(sc, PolygonRDDInputLocation, PolygonRDDStartOffset, PolygonRDDEndOffset, PolygonRDDSplitter, true);
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);

        objectRDD.spatialPartitioning(joinQueryPartitioningType);
        queryWindowRDD.spatialPartitioning(objectRDD.getPartitioner());

        objectRDD.buildIndex(PointRDDIndexType, true);

        objectRDD.indexedRDD.persist(StorageLevel.MEMORY_ONLY());
        queryWindowRDD.spatialPartitionedRDD.persist(StorageLevel.MEMORY_ONLY());

        for (int i = 0; i < eachQueryLoopTimes; i++) {
            long resultSize = JoinQuery.SpatialJoinQuery(objectRDD, queryWindowRDD, true, false).count();
            assert resultSize > 0;
        }
    }

    /**
     * Test spatial join query.
     *
     * @throws Exception the exception
     */
    public static void testDistanceJoinQuery()
            throws Exception
    {
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);
        CircleRDD queryWindowRDD = new CircleRDD(objectRDD, 0.1);

        objectRDD.spatialPartitioning(GridType.QUADTREE);
        queryWindowRDD.spatialPartitioning(objectRDD.getPartitioner());

        objectRDD.spatialPartitionedRDD.persist(StorageLevel.MEMORY_ONLY());
        queryWindowRDD.spatialPartitionedRDD.persist(StorageLevel.MEMORY_ONLY());

        for (int i = 0; i < eachQueryLoopTimes; i++) {

            long resultSize = JoinQuery.DistanceJoinQuery(objectRDD, queryWindowRDD, false, true).count();
            assert resultSize > 0;
        }
    }

    /**
     * Test spatial join query using index.
     *
     * @throws Exception the exception
     */
    public static void testDistanceJoinQueryUsingIndex()
            throws Exception
    {
        objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);
        CircleRDD queryWindowRDD = new CircleRDD(objectRDD, 0.1);

        objectRDD.spatialPartitioning(GridType.QUADTREE);
        queryWindowRDD.spatialPartitioning(objectRDD.getPartitioner());

        objectRDD.buildIndex(IndexType.RTREE, true);

        objectRDD.indexedRDD.persist(StorageLevel.MEMORY_ONLY());
        queryWindowRDD.spatialPartitionedRDD.persist(StorageLevel.MEMORY_ONLY());

        for (int i = 0; i < eachQueryLoopTimes; i++) {
            long resultSize = JoinQuery.DistanceJoinQuery(objectRDD, queryWindowRDD, true, true).count();
            assert resultSize > 0;
        }
    }

    public static void testLoadShapefileIntoPolygonRDD()
            throws Exception
    {
        ShapefileRDD shapefileRDD = new ShapefileRDD(sc, ShapeFileInputLocation);
        PolygonRDD spatialRDD = new PolygonRDD(shapefileRDD.getPolygonRDD());
        try {
            RangeQuery.SpatialRangeQuery(spatialRDD, new Envelope(-180, 180, -90, 90), false, false).count();
        }
        catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
}
