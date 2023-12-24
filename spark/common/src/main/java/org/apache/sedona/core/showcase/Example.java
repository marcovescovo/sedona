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
import org.locationtech.jts.geom.Polygon;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

// TODO: Auto-generated Javadoc

/**
 * The Class Example.
 */
public class Example
        implements Serializable
{

	public static BigInteger globalCounter = new BigInteger("0");
	public static int powFactor = 1;

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
     * The input location.
     */
    static String InputLocation;

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
     * The object RDD.
     */
    static PolygonRDD objectPolRDD;

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
     * The k NN query point.
     */
    static Polygon kNNQueryPolygon;

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
        SparkConf conf = new SparkConf().setAppName("SedonaRunnableExample").setMaster("local[1]");
        conf.set("spark.serializer", KryoSerializer.class.getName());
        conf.set("spark.kryo.registrator", SedonaKryoRegistrator.class.getName());

        sc = new JavaSparkContext(conf);
        Logger.getLogger("org").setLevel(Level.WARN);
        Logger.getLogger("akka").setLevel(Level.WARN);

        String resourceFolder = System.getProperty("user.dir") + "/src/test/resources/";

        PointRDDInputLocation = resourceFolder + "arealm-small.csv";
        InputLocation = resourceFolder + "primaryroads-polygon.csv";
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
        kNNQueryPoint = geometryFactory.createPoint(new Coordinate(0.5 * powFactor, 0.5 * powFactor));
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
    	
    	for (int i = 0; i <= 1000; i++) {
    		if (i != 1 && i != 10 && i != 50 && i != 100 && i != 1000) {
    			continue;
    		}
    		
    		for (int j = 1; j<=3; j++) {
    			
    			globalCounter = new BigInteger("0");
    	        
    	        try (BufferedReader br = new BufferedReader(new FileReader(InputLocation))) {

    	            String line;
    	            objectRDD = new PointRDD(sc, PointRDDInputLocation, PointRDDOffset, PointRDDSplitter, true);
    	            objectPolRDD = new PolygonRDD(sc, PolygonRDDInputLocation, PolygonRDDStartOffset, PolygonRDDNumPartitions, PolygonRDDSplitter, true);
    	            objectPolRDD.rawSpatialRDD.persist(StorageLevel.MEMORY_ONLY());
    	            objectPolRDD.buildIndex(PointRDDIndexType, false);
    	            LocalDateTime startDateTime = LocalDateTime.now();
    	            String formattedTimestamp = startDateTime.format(formatter);

    	            HashMap<Point,ArrayList> mappa = new HashMap<Point,ArrayList>();
    	            HashMap<Polygon,ArrayList> mappa_polygon = new HashMap<Polygon,ArrayList>();
    	            
    	            System.out.println("Test "+ i +"/"+ j +" - Start Timestamp: " + formattedTimestamp);
    	            
    	            while ((line = br.readLine()) != null) {
    	                ArrayList<Point> result_final = new ArrayList<Point>();
    	                ArrayList<Polygon> result_final_polygon = new ArrayList<Polygon>();
    	                // Process each line
    	                String[] fields = line.split(",");
    	                
    	                // Do something with the fields (for example, print them)
    	                for (String field : fields) {
    	                    kNNQueryPoint = geometryFactory.createPoint(new Coordinate(Double.parseDouble(fields[0]), Double.parseDouble(fields[1])));
    	                    kNNQueryPolygon = geometryFactory.createPolygon(geometryFactory.createLinearRing(new Coordinate[]{
    	                            new Coordinate(Double.parseDouble(fields[0]) * powFactor, Double.parseDouble(fields[1]) * powFactor),
    	                            new Coordinate(Double.parseDouble(fields[0])*2 * powFactor, Double.parseDouble(fields[1]) * powFactor),
    	                            new Coordinate(Double.parseDouble(fields[0])*2 * powFactor, Double.parseDouble(fields[1])*2 * powFactor),
    	                            new Coordinate(Double.parseDouble(fields[0]) * powFactor, Double.parseDouble(fields[1])*2 * powFactor),
    	                            new Coordinate(Double.parseDouble(fields[0]) * powFactor, Double.parseDouble(fields[1]) * powFactor)
    	                            }), null);
    	                	List<Polygon/*Point*/> result = KNNQuery.SpatialKnnQuery(objectPolRDD/*objectRDD*/, kNNQueryPolygon/*kNNQueryPoint*/, i, true);
    	                	//List<Polygon/*Point*/> result = KNNQuery.SpatialKnnQuery(objectPolRDD/*objectRDD*/, kNNQueryPolygon/*kNNQueryPoint*/, i, false);
    	                    assert result.size() > -1;
    	                    /*if (result.size() > 0) {
    	                    	for (int j = 0; j < result.size(); j++) {
    	                    		if (result_final_polygon.contains(result.get(j)) == false) {}
    	                            result_final_polygon.add(result.get(j));	
    	                    	}	
    	                    }
    	                    //result_final_polygon.sort(null);
    	                    if (result_final_polygon.size() > 0) {
    	                    	mappa_polygon.put(kNNQueryPolygon, result_final_polygon);
    	                    }*/
    	                    
    	                    //System.out.println(result_final.size());
    	                    //System.out.println();
    	                }
    	            }
    	        	LocalDateTime endDateTime = LocalDateTime.now();
    	            formattedTimestamp = endDateTime.format(formatter);
    	            System.out.println("End Timestamp: " + formattedTimestamp);
    	            System.out.println("Time elapsed: " + formatDuration(startDateTime, endDateTime));
    	            System.out.println("Operations performed: " + globalCounter.toString());
    	            System.out.println();

    	        } catch (IOException e) {
    	            e.printStackTrace();
    	        }
    			
    		}
    	}
        ////////////////////////////////////////////////////////////////
    }
    
    private static String formatDuration(LocalDateTime startTime, LocalDateTime endTime) {
    	
    	Duration duration = Duration.between(startTime, endTime);
    	
        long hours = duration.toHours();
        long minutes = duration.toMinutesPart();
        long seconds = duration.toSecondsPart();
        long millis = duration.toMillisPart();

        return String.format("%02d:%02d:%02d.%03d", hours, minutes, seconds, millis);
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
            List<Point> result = KNNQuery.SpatialKnnQuery(objectRDD, kNNQueryPoint, 10, true);
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
