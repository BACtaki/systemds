package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.lang.NotImplementedException;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.util.HashMap;

public class HyperLogLogHashBucket {
    private static final Log LOGGER = LogFactory.getLog(HyperLogLogHashBucket.class.getName());

    private static int N_BITS_HASH_BUCKET_KEY;
    private static HashMap<Integer, Integer> HASH_BUCKETS;
    private static final double correction = 0.77351;

    // TODO HLL
    // private TreeSet<Integer> sortedHashes;

    public HyperLogLogHashBucket(int leadingK) {
        N_BITS_HASH_BUCKET_KEY = leadingK;
        HASH_BUCKETS = new HashMap<>(N_BITS_HASH_BUCKET_KEY);
    }

    public void add(int hash) {
        // A | B
        // [1 1 ... 1] | [1 1 ... 1]
        // get N bits from left to create A and B
        int A = extractKBitsFromIndex(hash, N_BITS_HASH_BUCKET_KEY, 10);
        LOGGER.debug("A (decimal)=" + A + " A (binary)=" + Integer.toBinaryString(A));
        int B = extractKBitsFromIndex(hash, 10, 0);
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

        int m = HASH_BUCKETS.size();
        float sum = HASH_BUCKETS.values().stream().reduce(0, Integer::sum);

        return (int) Math.round((m * Math.pow(2, Math.floor(sum / m))) / correction);
    }

    public int getSuperLogLogEstimate() {
        throw new NotImplementedException("");
    }

    public int getHyperLogLogEstimate() {
        throw new NotImplementedException("");
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
