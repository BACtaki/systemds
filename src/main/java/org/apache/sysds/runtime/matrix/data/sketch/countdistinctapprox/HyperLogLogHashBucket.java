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

    private static int N_BITS_HASH_BUCKET_KEY;
    private static HashMap<Integer, Integer> HASH_BUCKETS;
    private static final double correction = 0.77351;

    public HyperLogLogHashBucket(int leadingK) {
        N_BITS_HASH_BUCKET_KEY = leadingK;
        HASH_BUCKETS = new HashMap<>(N_BITS_HASH_BUCKET_KEY);
    }

    public void add(int hash) {
        // A | B
        // [1 1 ... 1] | [1 1 ... 1]
        // get N bits from left to create A and B
        int A = extractKBitsFromIndex(hash, N_BITS_HASH_BUCKET_KEY, 32 - N_BITS_HASH_BUCKET_KEY);
        LOGGER.debug("A (decimal)=" + A + " A (binary)=" + Integer.toBinaryString(A));
        int B = extractKBitsFromIndex(hash, 32 - N_BITS_HASH_BUCKET_KEY, 0);
        LOGGER.debug("B (decimal)=" + B + " B (binary)=" + Integer.toBinaryString(B));

        // Use A to get bucket hash - x
        int existingZeroCount = HASH_BUCKETS.getOrDefault(A, 0);

        // Count number of leading 0s in B - y
        int newZeroCount = Integer.numberOfLeadingZeros(B);
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

        List<Integer> smallestKValues = getSmallestK(0.3);
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
     *
     * @param num
     * @param k
     * @param rightStartIndex 0-indexed from the right
     * @return
     */
    private int extractKBitsFromIndex(int num, int k, int rightStartIndex) {
        return ((1 << k) - 1) & (num >>> rightStartIndex);
    }
}
