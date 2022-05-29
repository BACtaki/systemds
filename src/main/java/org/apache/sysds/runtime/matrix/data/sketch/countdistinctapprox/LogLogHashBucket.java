package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.sysds.runtime.matrix.data.Pair;
import org.jetbrains.annotations.NotNull;

import java.util.*;

// Todo Parameterize?
public class LogLogHashBucket implements Iterable<Pair<Integer, Integer>> {

    public enum MeanType {
        SIMPLE, HARMONIC
    }

    private static final Log LOGGER = LogFactory.getLog(LogLogHashBucket.class.getName());

    // Number of bits used to map hash to bucket
    private static int K;
    // Maximum number of hash buckets
    private static int M;

    // Todo do we need to make this static? Pros and cons
    // Todo do we need to make this final? Pros and cons
    private HashMap<Integer, Integer> hashBuckets;
    private static final double correction = 0.77351;

    public LogLogHashBucket(int leadingK) {
        K = leadingK;
        M = (int) Math.pow(2, K);
        hashBuckets = new HashMap<>(M);
    }

    @NotNull
    @Override
    public Iterator<Pair<Integer, Integer>> iterator() {
        return new LogLogHashBucketIterator<>(hashBuckets);
    }

    public int size() {
        return hashBuckets.size();
    }

    /**
     * Follow method signature of Collection.add()
     * @param hash
     * @return
     */
    // Todo better method name
    public boolean addHash(int hash) {
        // A | B
        // [1 1 ... 1] | [1 1 ... 1]
        // get N bits from left to create A and B
        int A = extractKBitsFromRightIndex(hash, 32 - K, K);
        LOGGER.debug("A (decimal)=" + A + " A (binary)=" + Integer.toBinaryString(A));
        int B = extractKBitsFromRightIndex(hash, 0, 32 - K);
        LOGGER.debug("B (decimal)=" + B + " B (binary)=" + Integer.toBinaryString(B));

        // Use A to get bucket hash -> x
        int existingZeroCount = hashBuckets.getOrDefault(A, 0);

        // Position of leading 1 in binary representation of B -> y
        int newZeroCount = getBinaryRank(B);
        int longestSeqOfZeros = Math.max(existingZeroCount, newZeroCount);
        LOGGER.debug("Longest seq of 0s for key " + A + "=" + longestSeqOfZeros);

        // Now set BUCKETS[x] = y;
        return hashBuckets.put(A, longestSeqOfZeros) != null;
    }

    // Todo better method name
    public boolean put(int key, int value) {

        // Todo Replace value with maximum

        return hashBuckets.put(key, value) != null;
    }

    public void clear() {
        hashBuckets.clear();
    }

    public int getLogLogEstimate() {

        List<Integer> allValues = new LinkedList<>(hashBuckets.values());
        int m = allValues.size();

        int averageBitLengthFloor = (int)Math.floor(getMean(allValues, MeanType.SIMPLE));
        return (int) Math.round((m * Math.pow(2, averageBitLengthFloor)) / correction);
    }

    public int getSuperLogLogEstimate() {

        List<Integer> smallestKValues = getSmallestPercentage(0.7);
        int m = smallestKValues.size();

        int averageBitLengthFloor = (int)Math.floor(getMean(smallestKValues, MeanType.SIMPLE));
        return (int) Math.round((m * Math.pow(2, averageBitLengthFloor)) / correction);
    }

    public int getHyperLogLogEstimate() {

//        List<Integer> smallestKValues = getSmallestK(0.7);
        List<Integer> smallestKValues = getSmallestPercentage(0.7);
        int m = smallestKValues.size();

        int averageBitLengthFloor = (int)Math.floor(getMean(smallestKValues, MeanType.HARMONIC));
        return (int) Math.round((m * Math.pow(2, averageBitLengthFloor)) / correction);
    }

    public static double getMean(List<Integer> list, MeanType type) {
        int n = list.size();

        // Todo check if this makes sense
        if (n < 1) {
            return 0.0;
        }

        if (type == MeanType.HARMONIC) {
            double denominator = 0.0;
            for (Integer i : list) {
                // Todo division by zero error
                denominator += 1.0 / i;
            }
            return n / denominator;

        } else {  // default to simple mean
            return list.stream().reduce(0, Integer::sum) / (double)n;
        }

    }

    private List<Integer> getSmallestPercentage(double percentage) {
        LinkedList<Integer> values = new LinkedList<>(hashBuckets.values());
        Collections.sort(values);

        int removeN = (int) Math.floor((1 - percentage) * hashBuckets.size());
        for (int i=0; i<removeN; ++i) {
            values.pollLast();
        }

        return values;
    }

//    public void serialize(MatrixBlock blkOut, Types.Direction dir) {
//        int i = 0;
//        // HASH_BUCKETS will have a maximum of m = 2^k = 1024 entries
//        for (Integer hash : hashBuckets.keySet()) {
//            if (dir == Types.Direction.RowCol) {
//                // blkOut is a M x 2 matrix
//                blkOut.setValue(i, 0, hash);
//                blkOut.setValue(i, 1, hashBuckets.get(hash));
//            } else if (dir == Types.Direction.Row) {
//                // Todo fix
//                // blkOut is a M x 2 matrix
//                blkOut.setValue(i, 0, hash);
//                blkOut.setValue(i, 1, hashBuckets.get(hash));
//            } else {
//                // Todo fix
//                // blkOut is a 2 x N matrix
//                blkOut.setValue(0, i, hash);
//                blkOut.setValue(1, i, hashBuckets.get(hash));
//            }
//
//            i++;
//        }
//    }

    /**
     * Todo
     * The rank of a binary number is defined as the 1-indexed position of the leading 1 from the left?
     * @param hash
     * @return
     */
    private int getBinaryRank(int hash) {
        return Integer.numberOfLeadingZeros(hash) - K + 2;
    }

    /**
     *
     * @param num
     * @param rightIndex 0-indexed from the right
     * @param k number of bits
     * @return
     */
    private int extractKBitsFromRightIndex(int num, int rightIndex, int k) {
        return ((1 << k) - 1) & (num >>> rightIndex);
    }
}
