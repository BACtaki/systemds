package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

public class HyperLogLogHashBucket {
    private enum MeanType {
        SIMPLE, HARMONIC
    }

    private static final Log LOGGER = LogFactory.getLog(HyperLogLogHashBucket.class.getName());

    // Number of bits used to map hash to bucket
    private static int K;
    // Maximum number of hash buckets
    private static int M;

    private static HashMap<Integer, Integer> HASH_BUCKETS;
    private static final double correction = 0.77351;

    public HyperLogLogHashBucket(int leadingK) {
        K = leadingK;
        M = (int) Math.pow(2, K);
        HASH_BUCKETS = new HashMap<>(M);
    }

    public void add(int hash) {
        // A | B
        // [1 1 ... 1] | [1 1 ... 1]
        // get N bits from left to create A and B
        int A = extractKBitsFromRightIndex(hash, 32 - K, K);
        LOGGER.debug("A (decimal)=" + A + " A (binary)=" + Integer.toBinaryString(A));
        int B = extractKBitsFromRightIndex(hash, 0, 32 - K);
        LOGGER.debug("B (decimal)=" + B + " B (binary)=" + Integer.toBinaryString(B));

        // Use A to get bucket hash -> x
        int existingZeroCount = HASH_BUCKETS.getOrDefault(A, 0);

        // Position of leading 1 in binary representation of B -> y
        int newZeroCount = getBinaryRank(B);
        int longestSeqOfZeros = Math.max(existingZeroCount, newZeroCount);
        LOGGER.debug("Longest seq of 0s for key " + A + "=" + longestSeqOfZeros);

        // Now set BUCKETS[x] = y;
        HASH_BUCKETS.put(A, longestSeqOfZeros);
    }

    public int getLogLogEstimate() {

        List<Integer> allValues = new LinkedList<>(HASH_BUCKETS.values());
        int m = allValues.size();

        int averageBitLengthFloor = (int)Math.floor(getMean(allValues, MeanType.SIMPLE));
        return (int) Math.round((m * Math.pow(2, averageBitLengthFloor)) / correction);
    }

    public int getSuperLogLogEstimate() {

        List<Integer> smallestKValues = getSmallestK(0.7);
        int m = smallestKValues.size();

        int averageBitLengthFloor = (int)Math.floor(getMean(smallestKValues, MeanType.SIMPLE));
        return (int) Math.round((m * Math.pow(2, averageBitLengthFloor)) / correction);
    }

    public int getHyperLogLogEstimate() {

//        List<Integer> smallestKValues = getSmallestK(0.7);
        List<Integer> smallestKValues = getSmallestK(0.825);
        int m = smallestKValues.size();

        int averageBitLengthFloor = (int)Math.floor(getMean(smallestKValues, MeanType.HARMONIC));
        return (int) Math.round((m * Math.pow(2, averageBitLengthFloor)) / correction);
    }

    private double getMean(List<Integer> list, MeanType type) {
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

    private List<Integer> getSmallestK(double kPercentage) {
        LinkedList<Integer> values = new LinkedList<>(HASH_BUCKETS.values());
        Collections.sort(values);

        int removeN = (int) Math.floor((1 - kPercentage) * HASH_BUCKETS.size());
        for (int i=0; i<removeN; ++i) {
            values.pollLast();
        }

        return values;
    }

    /**
     * Todo
     * The rank of a binary number is defined as the 1-indexed position of the leading 1 from the left?
     * @param hash
     * @return
     */
    private int getBinaryRank(int hash) {
        return Integer.numberOfLeadingZeros(hash) - K + 1;
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
