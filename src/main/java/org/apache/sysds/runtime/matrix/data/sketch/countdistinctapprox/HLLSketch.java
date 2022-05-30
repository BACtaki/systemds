package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.lang.NotImplementedException;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.sysds.common.Types;
import org.apache.sysds.runtime.instructions.spark.data.CorrMatrixBlock;
import org.apache.sysds.runtime.matrix.data.MatrixBlock;
import org.apache.sysds.runtime.matrix.data.Pair;
import org.apache.sysds.runtime.matrix.operators.CountDistinctOperator;
import org.apache.sysds.utils.Hash;

public class HLLSketch extends CountDistinctApproxSketch {

    private static final Log LOGGER = LogFactory.getLog(HLLSketch.class.getName());

    // Todo rethink this comment
    // This is the number used to group the input values into "buckets".
    // Each input block will contain at most 1000 < 2^10 unique values, which means that
    // 10 bits are sufficient to count all unique numbers in any given matrix block.
    // The trailing bits are used to count, which leaves 32 - 10 = 22 bits to to used to determine
    // the bucket that the hash belongs to. Each bucket stores a 32bit int, so all buckets take a total of
    // 22 * 32 = 700 bits < 100 bytes of memory.

    // TODO Explore implications of using all 10 bits; do we even need to use all of them?
    // Choice of K:
    // In the CP case, blkIn is of size 1000 x 1000 -> there can be at most 10^6 unique values
    // (int) Math.ceil(Math.log(Math.pow(10, 6)) / Math.log(2)) = 20
    // Number of bits used to map hash to bucket
    private static final int K = 10;
    private static final double percentage = 0.7;
    private static final double correction = 0.77351;

    public HLLSketch(CountDistinctOperator op) {
        super(op);
    }

    @Override
    public Integer getScalarValue(MatrixBlock blkIn) {
        LogLogHashBucket hashBuckets = new LogLogHashBucket(K);
        for (int i=0; i<blkIn.getNumRows(); ++i) {
            for (int j=0; j<blkIn.getNumColumns(); ++j) {
                int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));

                hashBuckets.addHash(hash);
            }
        }

//        int estimate = hashBuckets.getLogLogEstimate();
//        int estimate = hashBuckets.getSuperLogLogEstimate();
        int estimate = hashBuckets.getHyperLogLogEstimate();

        if (blkIn.getNonZeros() != 0 && blkIn.getNonZeros() < estimate) {
            estimate = (int) blkIn.getNonZeros();
        }

        return estimate;
    }

    @Override
    public MatrixBlock getMatrixValue(CorrMatrixBlock arg0) {
        // Todo fix
        MatrixBlock blkIn = arg0.getValue();
        MatrixBlock blkInCorr = arg0.getCorrection();
        if (op.getDirection() == Types.Direction.RowCol) {
            int N = (int)blkInCorr.getValue(0, 0);

            LogLogHashBucket hashBuckets = new LogLogHashBucket(K);
            for (int i=0; i<N; ++i) {
                int k = (int) blkIn.getValue(i, 0);
                int v = (int) blkIn.getValue(i, 1);
                hashBuckets.put(k, v);
            }

            int estimate = hashBuckets.getHyperLogLogEstimate();

            // blkIn cannot be empty - reuse it
            MatrixBlock blkOut = blkIn.slice(0, 0, 0, 0);
            blkOut.setValue(0, 0, estimate);
            return blkOut;

        } else if (op.getDirection() == Types.Direction.Row) {

        } else {
            // Col
        }

        return new MatrixBlock(1);
    }

    @Override
    public CorrMatrixBlock create(MatrixBlock blkIn) {
        /**
         * blkIn - M x N, max 1000 x 1000
         *
         * The return value is a CorrMatrixBlock -> (HashBucket values, hash bucket metadata)
         *
         * Overwrite input matrix to store hash bucket values
         * N key-value pairs
         *
         * metadata ->
         *  n - number of hashes
         *  Todo parameterize this: smallestPercentage -> the smallest percentage of hashes to consider
         *  Currently, it is a constant that is applied at the very end (or should be).
         *
         */

        LogLogHashBucket hashBuckets = new LogLogHashBucket(K);
        int R = blkIn.getNumRows(), C = blkIn.getNumColumns();
        if (this.op.getDirection() == Types.Direction.RowCol) {
            int N = (int)Math.pow(2, K) + 1;

            MatrixBlock blkOutCorr = new MatrixBlock(2, N, false);
            // Push all values into the hash bucket data structure
            for (int i=0; i<R; ++i) {
                for (int j=0; j<C; ++j) {
                    int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                    LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));
                    hashBuckets.addHash(hash);
                }
            }

            // Serialize values in hash bucket data structure to a MatrixBlock
            int j = 0;
            for (Pair<Integer, Integer> kvPair : hashBuckets) {
                int k = kvPair.getKey();
                int v = kvPair.getValue();

                blkOutCorr.setValue(0, j, k);
                blkOutCorr.setValue(1, j, v);
                ++j;
            }

            blkOutCorr.setValue(0, N-1, hashBuckets.size());
            blkOutCorr.setValue(1, N-1, hashBuckets.size());

            return new CorrMatrixBlock(null, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Row) {
            int N = C + 1;

            MatrixBlock blkOut = blkIn;
            MatrixBlock blkOutCorr = new MatrixBlock(R, N, false);

            for (int i=0; i<R; ++i) {
                for (int j=0; j<C; ++j) {
                    int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                    LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));
                    hashBuckets.addHash(hash);
                }

                // Serialize hash buckets to MatrixBlock sketch
                int j = 0;
                for (Pair<Integer, Integer> kvPair : hashBuckets) {
                    int k = kvPair.getKey();
                    int v = kvPair.getValue();

                    blkOut.setValue(i, j, k);
                    blkOutCorr.setValue(i, j, v);
                    ++j;
                }
                blkOutCorr.setValue(i, N-1, hashBuckets.size());

                // Prepare hash buckets data structure for new row
                hashBuckets.clear();
            }

            return new CorrMatrixBlock(blkOut, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Col) {
            int N = R + 1;

            MatrixBlock blkOut = blkIn;
            MatrixBlock blkOutCorr = new MatrixBlock(N, C, false);

            for (int j=0; j<C; ++j) {
                for (int i=0; i<R; ++i) {
                    int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                    LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));
                    hashBuckets.addHash(hash);
                }

                // Serialize hash buckets to MatrixBlock sketch
                int i = 0;
                for (Pair<Integer, Integer> kvPair : hashBuckets) {
                    int k = kvPair.getKey();
                    int v = kvPair.getValue();

                    blkOut.setValue(i, j, k);
                    blkOutCorr.setValue(i, j, v);
                    ++i;
                }
                blkOutCorr.setValue(N-1, j, hashBuckets.size());

                // Prepare hash buckets data structure for new col
                hashBuckets.clear();
            }

            return new CorrMatrixBlock(blkOut, blkOutCorr);

        } else {
            throw new IllegalArgumentException("Invalid direction");
        }
    }

    @Override
    public CorrMatrixBlock union(CorrMatrixBlock arg0, CorrMatrixBlock arg1) {
        MatrixBlock blkInA = arg0.getValue();
        MatrixBlock blkInCorrA = arg0.getCorrection();

        MatrixBlock blkInB = arg1.getValue();
        MatrixBlock blkInCorrB = arg1.getCorrection();
        if (this.op.getDirection() == Types.Direction.RowCol) {
            // Both blkInCorrA and blkInCorrB are the same dimensions.
            // We will overwrite blkInCorrA to store the output.
            MatrixBlock blkOutCorr = blkInCorrA;

            LogLogHashBucket hashBuckets = new LogLogHashBucket(K);

            int nA = (int)blkInCorrA.getValue(0, blkInCorrA.getNumColumns() - 1);
            int nB = (int)blkInCorrB.getValue(0, blkInCorrB.getNumColumns() - 1);

            for (int j=0; j<Math.min(nA, nB); ++j) {
                int k = (int) blkInCorrA.getValue(0, j);
                int v = (int) blkInCorrA.getValue(1, j);
                hashBuckets.put(k, v);

                k = (int) blkInCorrB.getValue(0, j);
                v = (int) blkInCorrB.getValue(1, j);
                hashBuckets.put(k, v);
            }

            MatrixBlock larger = blkInA;
            if (Math.max(nA, nB) == nB) {
                larger = blkInB;
            }

            int j = Math.min(nA, nB);
            while (j < Math.max(nA, nB)) {
                int k = (int) larger.getValue(0, j);
                int v = (int) larger.getValue(1, j);
                hashBuckets.put(k, v);
                ++j;
            }

            j = 0;
            for (Pair<Integer, Integer> kvPair : hashBuckets) {
                int k = kvPair.getKey();
                int v = kvPair.getValue();

                blkOutCorr.setValue(0, j, k);
                blkOutCorr.setValue(1, j, k);
                ++j;
            }

            blkOutCorr.setValue(0, blkOutCorr.getNumColumns() - 1, hashBuckets.size());
            blkOutCorr.setValue(1, blkOutCorr.getNumColumns() - 1, hashBuckets.size());

            return new CorrMatrixBlock(null, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Row) {

            int nA = (int)blkInCorrA.getValue(0, blkInCorrA.getNumColumns() - 1);
            int nB = (int)blkInCorrB.getValue(0, blkInCorrB.getNumColumns() - 1);

            MatrixBlock blkOut = blkInA;
            MatrixBlock blkOutCorr = blkInCorrA;
            if (nB > nA) {
                blkOut = blkInB;
                blkOutCorr = blkInCorrB;
            }

            LogLogHashBucket hashBuckets = new LogLogHashBucket(K);

            return new CorrMatrixBlock(null, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Col) {

        } else {
            throw new IllegalArgumentException("Invalid direction");
        }

        throw new NotImplementedException("union has not been implemented yet");
    }

    @Override
    public CorrMatrixBlock intersection(CorrMatrixBlock arg0, CorrMatrixBlock arg1) {
        throw new NotImplementedException(String.format("%s intersection has not been implemented yet",
                HLLSketch.class.getSimpleName()));
    }
}
