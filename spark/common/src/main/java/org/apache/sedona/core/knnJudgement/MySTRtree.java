/*
 * Copyright (c) 2016 Vivid Solutions.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */
package org.apache.sedona.core.knnJudgement;

import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

import org.apache.sedona.core.showcase.Example;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.ItemVisitor;
import org.locationtech.jts.index.SpatialIndex;
import org.locationtech.jts.index.strtree.AbstractNode;
import org.locationtech.jts.index.strtree.AbstractSTRtree;
import org.locationtech.jts.index.strtree.Boundable;
import org.locationtech.jts.index.strtree.EnvelopeDistance;
import org.locationtech.jts.index.strtree.ItemBoundable;
import org.locationtech.jts.index.strtree.ItemDistance;
import org.locationtech.jts.util.Assert;

class BoundablePair implements Comparable {
	private Boundable boundable1;
	private Boundable boundable2;
	private double distance;
	private ItemDistance itemDistance;
	// private double maxDistance = -1.0;

	public BoundablePair(Boundable boundable1, Boundable boundable2, ItemDistance itemDistance) {
		this.boundable1 = boundable1;
		this.boundable2 = boundable2;
		this.itemDistance = itemDistance;
		distance = distance();
	}

	/**
	 * Gets one of the member {@link Boundable}s in the pair (indexed by [0, 1]).
	 * 
	 * @param i the index of the member to return (0 or 1)
	 * @return the chosen member
	 */
	public Boundable getBoundable(int i) {
		if (i == 0)
			return boundable1;
		return boundable2;
	}

	/**
	 * Computes the maximum distance between any two items in the pair of nodes.
	 * 
	 * @return the maximum distance between items in the pair
	 */
	public double maximumDistance() {
		return EnvelopeDistance.maximumDistance((Envelope) boundable1.getBounds(), (Envelope) boundable2.getBounds());
	}

	/**
	 * Computes the distance between the {@link Boundable}s in this pair. The
	 * boundables are either composites or leaves. If either is composite, the
	 * distance is computed as the minimum distance between the bounds. If both are
	 * leaves, the distance is computed by
	 * {@link #itemDistance(ItemBoundable, ItemBoundable)}.
	 * 
	 * @return
	 */
	private double distance() {
		// if items, compute exact distance
		if (isLeaves()) {
			return itemDistance.distance((ItemBoundable) boundable1, (ItemBoundable) boundable2);
		}
		// otherwise compute distance between bounds of boundables
		return ((Envelope) boundable1.getBounds()).distance(((Envelope) boundable2.getBounds()));
	}

	/**
	 * Gets the minimum possible distance between the Boundables in this pair. If
	 * the members are both items, this will be the exact distance between them.
	 * Otherwise, this distance will be a lower bound on the distances between the
	 * items in the members.
	 * 
	 * @return the exact or lower bound distance for this pair
	 */
	public double getDistance() {
		return distance;
	}

	/**
	 * Compares two pairs based on their minimum distances
	 */
	public int compareTo(Object o) {
		BoundablePair nd = (BoundablePair) o;
		if (distance < nd.distance)
			return -1;
		if (distance > nd.distance)
			return 1;
		return 0;
	}

	/**
	 * Tests if both elements of the pair are leaf nodes
	 * 
	 * @return true if both pair elements are leaf nodes
	 */
	public boolean isLeaves() {
		return !(isComposite(boundable1) || isComposite(boundable2));
	}

	public static boolean isComposite(Object item) {
		return (item instanceof AbstractNode);
	}

	private static double area(Boundable b) {
		return ((Envelope) b.getBounds()).getArea();
	}

	/**
	 * For a pair which is not a leaf (i.e. has at least one composite boundable)
	 * computes a list of new pairs from the expansion of the larger boundable with
	 * distance less than minDistance and adds them to a priority queue.
	 * <p>
	 * Note that expanded pairs may contain the same item/node on both sides. This
	 * must be allowed to support distance functions which have non-zero distances
	 * between the item and itself (non-zero reflexive distance).
	 * 
	 * @param priQ        the priority queue to add the new pairs to
	 * @param minDistance the limit on the distance between added pairs
	 * 
	 */
	public void expandToQueue(PriorityQueue priQ, double minDistance) {
		boolean isComp1 = isComposite(boundable1);
		boolean isComp2 = isComposite(boundable2);

		/**
		 * HEURISTIC: If both boundable are composite, choose the one with largest area
		 * to expand. Otherwise, simply expand whichever is composite.
		 */
		if (isComp1 && isComp2) {
			if (area(boundable1) > area(boundable2)) {
				expand(boundable1, boundable2, false, priQ, minDistance);
				return;
			} else {
				expand(boundable2, boundable1, true, priQ, minDistance);
				return;
			}
		} else if (isComp1) {
			expand(boundable1, boundable2, false, priQ, minDistance);
			return;
		} else if (isComp2) {
			expand(boundable2, boundable1, true, priQ, minDistance);
			return;
		}

		throw new IllegalArgumentException("neither boundable is composite");
	}

	private void expand(Boundable bndComposite, Boundable bndOther, boolean isFlipped, PriorityQueue priQ,
			double minDistance) {
		List children = ((AbstractNode) bndComposite).getChildBoundables();
		for (Iterator i = children.iterator(); i.hasNext();) {
			Boundable child = (Boundable) i.next();
			BoundablePair bp;
			if (isFlipped) {
				bp = new BoundablePair(bndOther, child, itemDistance);
			} else {
				bp = new BoundablePair(child, bndOther, itemDistance);
			}
			// only add to queue if this pair might contain the closest points
			// MD - it's actually faster to construct the object rather than called
			// distance(child, bndOther)!
			if (bp.getDistance() < minDistance) {
				priQ.add(bp);
			}
		}
	}
}

/**
 * A query-only R-tree created using the Sort-Tile-Recursive (STR) algorithm.
 * For two-dimensional spatial data.
 * <P>
 * The STR packed R-tree is simple to implement and maximizes space utilization;
 * that is, as many leaves as possible are filled to capacity. Overlap between
 * nodes is far less than in a basic R-tree. However, the index is semi-static;
 * once the tree has been built (which happens automatically upon the first
 * query), items may not be added. Items may be removed from the tree using
 * {@link #remove(Envelope, Object)}.
 * <P>
 * Described in: P. Rigaux, Michel Scholl and Agnes Voisard. <i>Spatial
 * Databases With Application To GIS</i>. Morgan Kaufmann, San Francisco, 2002.
 * <p>
 * <b>Note that inserting items into a tree is not thread-safe.</b> Inserting
 * performed on more than one thread must be synchronized externally.
 * <p>
 * Querying a tree is thread-safe. The building phase is done synchronously, and
 * querying is stateless.
 *
 * @version 1.7
 */
public class MySTRtree extends AbstractSTRtree implements SpatialIndex, Serializable {

	static final class STRtreeNode extends AbstractNode {
		STRtreeNode(int level) {
			super(level);
		}

		protected Object computeBounds() {
			Envelope bounds = null;
			for (Iterator i = getChildBoundables().iterator(); i.hasNext();) {
				Boundable childBoundable = (Boundable) i.next();
				if (bounds == null) {
					bounds = new Envelope((Envelope) childBoundable.getBounds());
				} else {
					bounds.expandToInclude((Envelope) childBoundable.getBounds());
				}
			}
			return bounds;
		}
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 259274702368956900L;

	private static Comparator xComparator = new Comparator() {
		public int compare(Object o1, Object o2) {
			return compareDoubles(centreX((Envelope) ((Boundable) o1).getBounds()),
					centreX((Envelope) ((Boundable) o2).getBounds()));
		}
	};
	private static Comparator yComparator = new Comparator() {
		public int compare(Object o1, Object o2) {
			return compareDoubles(centreY((Envelope) ((Boundable) o1).getBounds()),
					centreY((Envelope) ((Boundable) o2).getBounds()));
		}
	};

	private static double centreX(Envelope e) {
		return avg(e.getMinX(), e.getMaxX());
	}

	private static double centreY(Envelope e) {
		return avg(e.getMinY(), e.getMaxY());
	}

	private static double avg(double a, double b) {
		return (a + b) / 2d;
	}

	private static IntersectsOp intersectsOp = new IntersectsOp() {
		public boolean intersects(Object aBounds, Object bBounds) {
			return ((Envelope) aBounds).intersects((Envelope) bBounds);
		}
	};

	/**
	 * Creates the parent level for the given child level. First, orders the items
	 * by the x-values of the midpoints, and groups them into vertical slices. For
	 * each slice, orders the items by the y-values of the midpoints, and group them
	 * into runs of size M (the node capacity). For each run, creates a new (parent)
	 * node.
	 */
	protected List createParentBoundables(List childBoundables, int newLevel) {
		Assert.isTrue(!childBoundables.isEmpty());
		int minLeafCount = (int) Math.ceil((childBoundables.size() / (double) getNodeCapacity()));
		ArrayList sortedChildBoundables = new ArrayList(childBoundables);
		Collections.sort(sortedChildBoundables, xComparator);
		List[] verticalSlices = verticalSlices(sortedChildBoundables, (int) Math.ceil(Math.sqrt(minLeafCount)));
		return createParentBoundablesFromVerticalSlices(verticalSlices, newLevel);
	}

	private List createParentBoundablesFromVerticalSlices(List[] verticalSlices, int newLevel) {
		Assert.isTrue(verticalSlices.length > 0);
		List parentBoundables = new ArrayList();
		for (int i = 0; i < verticalSlices.length; i++) {
			parentBoundables.addAll(createParentBoundablesFromVerticalSlice(verticalSlices[i], newLevel));
		}
		return parentBoundables;
	}

	protected List createParentBoundablesFromVerticalSlice(List childBoundables, int newLevel) {
		return super.createParentBoundables(childBoundables, newLevel);
	}

	/**
	 * @param childBoundables Must be sorted by the x-value of the envelope
	 *                        midpoints
	 */
	protected List[] verticalSlices(List childBoundables, int sliceCount) {
		int sliceCapacity = (int) Math.ceil(childBoundables.size() / (double) sliceCount);
		List[] slices = new List[sliceCount];
		Iterator i = childBoundables.iterator();
		for (int j = 0; j < sliceCount; j++) {
			slices[j] = new ArrayList();
			int boundablesAddedToSlice = 0;
			while (i.hasNext() && boundablesAddedToSlice < sliceCapacity) {
				Boundable childBoundable = (Boundable) i.next();
				slices[j].add(childBoundable);
				boundablesAddedToSlice++;
			}
		}
		return slices;
	}

	private static final int DEFAULT_NODE_CAPACITY = 10;

	/**
	 * Constructs an STRtree with the default node capacity.
	 */
	public MySTRtree() {
		this(DEFAULT_NODE_CAPACITY);
	}

	/**
	 * Constructs an STRtree with the given maximum number of child nodes that a
	 * node may have.
	 * <p>
	 * The minimum recommended capacity setting is 4.
	 * 
	 */
	public MySTRtree(int nodeCapacity) {
		super(nodeCapacity);
	}

	/**
	 * Constructs an STRtree with the given maximum number of child nodes that a
	 * node may have, and the root that links to all other nodes
	 * <p>
	 * The minimum recommended capacity setting is 4.
	 *
	 */
	public MySTRtree(int nodeCapacity, STRtreeNode root) {
		super(nodeCapacity, root);
	}

	/**
	 * Constructs an STRtree with the given maximum number of child nodes that a
	 * node may have, and all leaf nodes in the tree
	 * <p>
	 * The minimum recommended capacity setting is 4.
	 *
	 */
	public MySTRtree(int nodeCapacity, ArrayList itemBoundables) {
		super(nodeCapacity, itemBoundables);
	}

	protected AbstractNode createNode(int level) {
		return new STRtreeNode(level);
	}

	protected IntersectsOp getIntersectsOp() {
		return intersectsOp;
	}

	/**
	 * Inserts an item having the given bounds into the tree.
	 */
	public void insert(Envelope itemEnv, Object item) {
		if (itemEnv.isNull()) {
			return;
		}
		super.insert(itemEnv, item);
	}

	/**
	 * Returns items whose bounds intersect the given envelope.
	 */
	public List query(Envelope searchEnv) {
		// Yes this method does something. It specifies that the bounds is an
		// Envelope. super.query takes an Object, not an Envelope. [Jon Aquino
		// 10/24/2003]
		return super.query((Object) searchEnv);
	}

	/**
	 * Returns items whose bounds intersect the given envelope.
	 */
	public void query(Envelope searchEnv, ItemVisitor visitor) {
		// Yes this method does something. It specifies that the bounds is an
		// Envelope. super.query takes an Object, not an Envelope. [Jon Aquino
		// 10/24/2003]
		super.query(searchEnv, visitor);
	}

	/**
	 * Removes a single item from the tree.
	 *
	 * @param itemEnv the Envelope of the item to remove
	 * @param item    the item to remove
	 * @return <code>true</code> if the item was found
	 */
	public boolean remove(Envelope itemEnv, Object item) {
		return super.remove(itemEnv, item);
	}

	/**
	 * Returns the number of items in the tree.
	 *
	 * @return the number of items in the tree
	 */
	public int size() {
		return super.size();
	}

	/**
	 * Returns the number of items in the tree.
	 *
	 * @return the number of items in the tree
	 */
	public int depth() {
		return super.depth();
	}

	protected Comparator getComparator() {
		return yComparator;
	}

	/**
	 * Finds the two nearest items in the tree, using {@link ItemDistance} as the
	 * distance metric. A Branch-and-Bound tree traversal algorithm is used to
	 * provide an efficient search.
	 * <p>
	 * If the tree is empty, the return value is <code>null</code.
	 * If the tree contains only one item, 
	 * the return value is a pair containing that item.  
	 * <b>
	 * If it is required to find only pairs of distinct items,
	 * the {@link ItemDistance} function must be <b>anti-reflexive</b>.
	 * 
	 * &#64;param itemDist a distance metric applicable to the items in this tree
	 * @return the pair of the nearest items
	 *    or <code>null</code> if the tree is empty
	 */
	public Object[] nearestNeighbour(ItemDistance itemDist) {
		if (isEmpty())
			return null;

		// if tree has only one item this will return null
		BoundablePair bp = new BoundablePair(this.getRoot(), this.getRoot(), itemDist);
		return nearestNeighbour(bp);
	}

	/**
	 * Finds the item in this tree which is nearest to the given {@link Object},
	 * using {@link ItemDistance} as the distance metric. A Branch-and-Bound tree
	 * traversal algorithm is used to provide an efficient search.
	 * <p>
	 * The query <tt>object</tt> does <b>not</b> have to be contained in the tree,
	 * but it does have to be compatible with the <tt>itemDist</tt> distance metric.
	 * 
	 * @param env      the envelope of the query item
	 * @param item     the item to find the nearest neighbour of
	 * @param itemDist a distance metric applicable to the items in this tree and
	 *                 the query item
	 * @return the nearest item in this tree or <code>null</code> if the tree is
	 *         empty
	 */
	public Object nearestNeighbour(Envelope env, Object item, ItemDistance itemDist) {
		if (isEmpty())
			return null;

		Boundable bnd = new ItemBoundable(env, item);
		BoundablePair bp = new BoundablePair(this.getRoot(), bnd, itemDist);
		return nearestNeighbour(bp)[0];
	}

	/**
	 * Finds the two nearest items from this tree and another tree, using
	 * {@link ItemDistance} as the distance metric. A Branch-and-Bound tree
	 * traversal algorithm is used to provide an efficient search. The result value
	 * is a pair of items, the first from this tree and the second from the argument
	 * tree.
	 * 
	 * @param tree     another tree
	 * @param itemDist a distance metric applicable to the items in the trees
	 * @return the pair of the nearest items, one from each tree or
	 *         <code>null</code> if no pair of distinct items can be found
	 */
	public Object[] nearestNeighbour(MySTRtree tree, ItemDistance itemDist) {
		if (isEmpty() || tree.isEmpty())
			return null;
		BoundablePair bp = new BoundablePair(this.getRoot(), tree.getRoot(), itemDist);
		return nearestNeighbour(bp);
	}

	private Object[] nearestNeighbour(BoundablePair initBndPair) {
		double distanceLowerBound = Double.POSITIVE_INFINITY;
		BoundablePair minPair = null;

		// initialize search queue
		PriorityQueue priQ = new PriorityQueue();
		priQ.add(initBndPair);

		while (!priQ.isEmpty() && distanceLowerBound > 0.0) {
			// pop head of queue and expand one side of pair
			BoundablePair bndPair = (BoundablePair) priQ.poll();
			double pairDistance = bndPair.getDistance();

			/**
			 * If the distance for the first pair in the queue is >= current minimum
			 * distance, other nodes in the queue must also have a greater distance. So the
			 * current minDistance must be the true minimum, and we are done.
			 */
			if (pairDistance >= distanceLowerBound)
				break;

			/**
			 * If the pair members are leaves then their distance is the exact lower bound.
			 * Update the distanceLowerBound to reflect this (which must be smaller, due to
			 * the test immediately prior to this).
			 */
			if (bndPair.isLeaves()) {
				// assert: currentDistance < minimumDistanceFound
				distanceLowerBound = pairDistance;
				minPair = bndPair;
			} else {
				/**
				 * Otherwise, expand one side of the pair, and insert the expanded pairs into
				 * the queue. The choice of which side to expand is determined heuristically.
				 */
				bndPair.expandToQueue(priQ, distanceLowerBound);
			}
		}
		if (minPair == null)
			return null;
		// done - return items with min distance
		return new Object[] { ((ItemBoundable) minPair.getBoundable(0)).getItem(),
				((ItemBoundable) minPair.getBoundable(1)).getItem() };
	}

	/**
	 * Tests whether some two items from this tree and another tree lie within a
	 * given distance. {@link ItemDistance} is used as the distance metric. A
	 * Branch-and-Bound tree traversal algorithm is used to provide an efficient
	 * search.
	 * 
	 * @param tree        another tree
	 * @param itemDist    a distance metric applicable to the items in the trees
	 * @param maxDistance the distance limit for the search
	 * @return true if there are items within the distance
	 */
	public boolean isWithinDistance(MySTRtree tree, ItemDistance itemDist, double maxDistance) {
		BoundablePair bp = new BoundablePair(this.getRoot(), tree.getRoot(), itemDist);
		return isWithinDistance(bp, maxDistance);
	}

	/**
	 * Performs a withinDistance search on the tree node pairs. This is a different
	 * search algorithm to nearest neighbour. It can utilize the
	 * {@link BoundablePair#maximumDistance()} between tree nodes to confirm if two
	 * internal nodes must have items closer than the maxDistance, and short-circuit
	 * the search.
	 * 
	 * @param initBndPair the initial pair containing the tree root nodes
	 * @param maxDistance the maximum distance to search for
	 * @return true if two items lie within the given distance
	 */
	private boolean isWithinDistance(BoundablePair initBndPair, double maxDistance) {
		double distanceUpperBound = Double.POSITIVE_INFINITY;

		// initialize search queue
		PriorityQueue priQ = new PriorityQueue();
		priQ.add(initBndPair);

		while (!priQ.isEmpty()) {
			// pop head of queue and expand one side of pair
			BoundablePair bndPair = (BoundablePair) priQ.poll();
			double pairDistance = bndPair.getDistance();

			/**
			 * If the distance for the first pair in the queue is > maxDistance, all other
			 * pairs in the queue must have a greater distance as well. So can conclude no
			 * items are within the distance and terminate with result = false
			 */
			if (pairDistance > maxDistance)
				return false;

			/**
			 * If the maximum distance between the nodes is less than the maxDistance, than
			 * all items in the nodes must be closer than the max distance. Then can
			 * terminate with result = true.
			 * 
			 * NOTE: using Envelope MinMaxDistance would provide a tighter bound, but not
			 * much performance improvement has been observed
			 */
			if (bndPair.maximumDistance() <= maxDistance)
				return true;
			/**
			 * If the pair items are leaves then their actual distance is an upper bound.
			 * Update the distanceUpperBound to reflect this
			 */
			if (bndPair.isLeaves()) {
				// assert: currentDistance < minimumDistanceFound
				distanceUpperBound = pairDistance;

				/**
				 * If the items are closer than maxDistance can terminate with result = true.
				 */
				if (distanceUpperBound <= maxDistance)
					return true;
			} else {
				/**
				 * Otherwise, expand one side of the pair, and insert the expanded pairs into
				 * the queue. The choice of which side to expand is determined heuristically.
				 */
				bndPair.expandToQueue(priQ, distanceUpperBound);
			}
		}
		return false;
	}

	/**
	 * Finds up to k items in this tree which are the nearest neighbors to the given
	 * {@code item}, using {@code itemDist} as the distance metric. A
	 * Branch-and-Bound tree traversal algorithm is used to provide an efficient
	 * search. This method implements the KNN algorithm described in the following
	 * paper:
	 * <p>
	 * Roussopoulos, Nick, Stephen Kelley, and Frédéric Vincent. "Nearest neighbor
	 * queries." ACM sigmod record. Vol. 24. No. 2. ACM, 1995.
	 * <p>
	 * The query {@code item} does <b>not</b> have to be contained in the tree, but
	 * it does have to be compatible with the {@code itemDist} distance metric.
	 * <p>
	 * If the tree size is smaller than k fewer items will be returned. If the tree
	 * is empty an array of size 0 is returned.
	 * 
	 * @param env      the envelope of the query item
	 * @param item     the item to find the nearest neighbours of
	 * @param itemDist a distance metric applicable to the items in this tree and
	 *                 the query item
	 * @param k        the maximum number of nearest items to search for
	 * @return an array of the nearest items found (with length between 0 and K)
	 */
	public Object[] nearestNeighbour(Envelope env, Object item, ItemDistance itemDist, int k, double d) {
		if (isEmpty())
			return new Object[0];

		Boundable bnd = new ItemBoundable(env, item);
		BoundablePair bp = new BoundablePair(this.getRoot(), bnd, itemDist);
		return nearestNeighbourK(bp, k, d);
	}

	private Object[] nearestNeighbourK(BoundablePair initBndPair, int k, double d) {
		return nearestNeighbourK(initBndPair, d, k);
	}

	private Object[] nearestNeighbourK(BoundablePair initBndPair, double maxDistance, int k) {
		
		double distanceLowerBound = maxDistance;

		// initialize internal structures
		PriorityQueue priQ = new PriorityQueue();

		// initialize queue
		priQ.add(initBndPair);

		PriorityQueue kNearestNeighbors = new PriorityQueue();

		BigInteger counter = BigInteger.ZERO;

		while (!priQ.isEmpty() && distanceLowerBound >= 0.0000) {

			counter = counter.add(new BigInteger("1"));

			// pop head of queue and expand one side of pair
			BoundablePair bndPair = (BoundablePair) priQ.poll();
			double pairDistance = bndPair.getDistance();

			/**
			 * If the distance for the first node in the queue is >= the current maximum
			 * distance in the k queue , all other nodes in the queue must also have a
			 * greater distance. So the current minDistance must be the true minimum, and we
			 * are done.
			 */
			if (pairDistance >= distanceLowerBound) {
				break;
			}
			/**
			 * If the pair members are leaves then their distance is the exact lower bound.
			 * Update the distanceLowerBound to reflect this (which must be smaller, due to
			 * the test immediately prior to this).
			 */
			if (bndPair.isLeaves()) {
				// assert: currentDistance < minimumDistanceFound

				if (kNearestNeighbors.size() < k) {
					kNearestNeighbors.add(bndPair);
				} else {

					BoundablePair bp1 = (BoundablePair) kNearestNeighbors.peek();
					if (bp1.getDistance() > pairDistance) {
						kNearestNeighbors.poll();
						kNearestNeighbors.add(bndPair);
					}
					/*
					 * minDistance should be the farthest point in the K nearest neighbor queue.
					 */
					BoundablePair bp2 = (BoundablePair) kNearestNeighbors.peek();
					distanceLowerBound = bp2.getDistance();
				}
			} else {
				/**
				 * Otherwise, expand one side of the pair, (the choice of which side to expand
				 * is heuristically determined) and insert the new expanded pairs into the queue
				 */
				bndPair.expandToQueue(priQ, distanceLowerBound);
			}
		}
		// done - return items with min distance
        Example.globalCounter = Example.globalCounter.add(counter);

		return getItems(kNearestNeighbors);
	}

	private static Object[] getItems(PriorityQueue kNearestNeighbors) {
		/**
		 * Iterate the K Nearest Neighbour Queue and retrieve the item from each
		 * BoundablePair in this queue
		 */
		Object[] items = new Object[kNearestNeighbors.size()];
		int count = 0;
		while (!kNearestNeighbors.isEmpty()) {
			BoundablePair bp = (BoundablePair) kNearestNeighbors.poll();
			items[count] = ((ItemBoundable) bp.getBoundable(0)).getItem();
			count++;
		}
		return items;
	}
}
