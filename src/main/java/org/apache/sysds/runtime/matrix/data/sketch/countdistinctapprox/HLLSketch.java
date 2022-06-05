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
        int R = blkIn.getNumRows();
        int C = blkIn.getNumColumns();

        if (this.op.getDirection() == Types.Direction.RowCol) {
            // Todo We need a better name for N
            int maxNumHashes = (int)Math.pow(2, K) + 1;

            MatrixBlock blkOutCorr = new MatrixBlock(2, maxNumHashes, false);
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

            // The exact number of hashes is stored in the last column of
            // the correction matrix in the RowCol case
            blkOutCorr.setValue(0, maxNumHashes-1, hashBuckets.size());
            blkOutCorr.setValue(1, maxNumHashes-1, hashBuckets.size());

            return new CorrMatrixBlock(null, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Row) {
            int maxNumHashes = C + 1;

            MatrixBlock blkOut = blkIn;
            MatrixBlock blkOutCorr = new MatrixBlock(R, maxNumHashes, false);

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

                // The exact number of hashes is stored in the last column of
                // each row in the correction matrix
                blkOutCorr.setValue(i, maxNumHashes-1, hashBuckets.size());

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

                // The exact number of hashes is stored in the last row of
                // each column in the correction matrix
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

        // Todo we have to be careful when taking the union in the row/col cases:
        //   we have to store the union in the larger of the 2 and ensure that
        //   reduceByKey() always retains the larger of the 2
        //   -> will swapping at the beginning of each case work?
        //   Swap to ensure the left matrix is always the larger of the 2
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
                blkOutCorr.setValue(1, j, v);
                ++j;
            }

            blkOutCorr.setValue(0, blkOutCorr.getNumColumns() - 1, hashBuckets.size());
            blkOutCorr.setValue(1, blkOutCorr.getNumColumns() - 1, hashBuckets.size());

            return new CorrMatrixBlock(null, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Row) {

            // Todo Proposed optimization: unroll loop (?) and iterate over 2 sketches simultaneously
            // [P] | [Q]
            // [R] | [S]
            // Assume [P] is overlapping sub-matrix.
            // The proposed optimization strategy:
            // 1. Iterate through both sketches simultaneously for [P]
            // 2. Iterate through larger sketch only for [Q], [R], and [S]
            //
            // It may be easier to iterate through [Q] + [S] and [R] + [S].
            // Double iteration over [S] will not affect correctness as the hash bucket
            // data structure retains the largest key anyway.
            //
            // It would affect performance though:
            // The difference in size between [P] and [S] will affect performance most because
            // the main objective of this optimization is to prevent iterating over 2 sub-matrices
            // of size [P] twice in series for the input matrices separately.

            LogLogHashBucket hashBuckets = new LogLogHashBucket(K);

            MatrixBlock tallerMatrix = blkInA;
            MatrixBlock tallerMatrixCorr = blkInCorrA;
            int largerR = blkInA.getNumRows();
            int smallerR = blkInB.getNumRows();
            if (blkInA.getNumRows() < blkInB.getNumRows()) {
                tallerMatrix = blkInB;
                tallerMatrixCorr = blkInCorrB;
                largerR = blkInB.getNumRows();
                smallerR = blkInA.getNumRows();
            }

            // Allocate new memory: Todo rationale
            MatrixBlock blkOut = new MatrixBlock(Math.max(blkInA.getNumRows(), blkInB.getNumRows()), Math.max(blkInA.getNumColumns(), blkInB.getNumColumns()), false);
            MatrixBlock blkOutCorr = new MatrixBlock(Math.max(blkInA.getNumRows(), blkInB.getNumRows()), Math.max(blkInA.getNumColumns(), blkInB.getNumColumns()) + 1, false);

            for (int i=0; i<smallerR; ++i) {
                int nA = (int)blkInCorrA.getValue(i, blkInCorrA.getNumColumns() - 1);
                int nB = (int)blkInCorrB.getValue(i, blkInCorrB.getNumColumns() - 1);

                MatrixBlock narrowerMatrix = blkInA;
                MatrixBlock narrowerMatrixCorr = blkInCorrA;
                int smallerC = blkInA.getNumColumns();
                int largerC = blkInB.getNumColumns();
                if (nB < nA) {
                    narrowerMatrix = blkInB;
                    narrowerMatrixCorr = blkInCorrB;
                    smallerC = blkInB.getNumColumns();
                    largerC = blkInA.getNumColumns();
                }

                for (int j=0; j<smallerC; ++j) {
                    int k = (int)blkInA.getValue(i, j);
                    int v = (int)blkInCorrA.getValue(i, j);
                    hashBuckets.put(k, v);

                    k = (int)blkInB.getValue(i, j);
                    v = (int)blkInCorrB.getValue(i, j);
                    hashBuckets.put(k, v);
                }

                // Todo fix bug: need wider matrix here instead of narrower
                for (int j=smallerC; j<largerC; ++j) {
                    int k = (int)narrowerMatrix.getValue(i, j);
                    int v = (int)narrowerMatrixCorr.getValue(i, j);
                    hashBuckets.put(k, v);
                }

                // Serialize hash bucket for this row to output matrix
                int j=0;
                for (Pair<Integer, Integer> kvPair : hashBuckets) {
                    int k = kvPair.getKey();
                    int v = kvPair.getValue();
                    blkOut.setValue(i, j, k);
                    blkOutCorr.setValue(i, j, v);
                }

                blkOutCorr.setValue(i, blkOutCorr.getNumColumns() - 1, hashBuckets.size());

                hashBuckets.clear();
            }

            // Copy out the taller matrix verbatim - no need for comparison
            for (int i=smallerR; i<largerR; ++i) {
                int nHashes = (int)tallerMatrix.getValue(i, tallerMatrix.getNumColumns() - 1);
                for (int j=0; j<nHashes; ++j) {
                    int k = (int)tallerMatrix.getValue(i, j);
                    int v = (int)tallerMatrixCorr.getValue(i, j);

                    blkOut.setValue(i, j, k);
                    blkOutCorr.setValue(i, j, v);
                }

                blkInCorrA.setValue(i, blkInCorrA.getNumColumns() - 1, nHashes);
            }

            return new CorrMatrixBlock(blkOut, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Col) {

            LogLogHashBucket hashBuckets = new LogLogHashBucket(K);

            MatrixBlock widerMatrix = blkInA;
            MatrixBlock widerMatrixCorr = blkInCorrA;
            int largerC = blkInA.getNumColumns();
            int smallerC = blkInB.getNumColumns();
            if (blkInA.getNumColumns() < blkInB.getNumColumns()) {
                widerMatrix = blkInB;
                widerMatrixCorr = blkInCorrB;
                largerC = blkInB.getNumColumns();
                smallerC = blkInA.getNumColumns();
            }

            MatrixBlock blkOut = new MatrixBlock(Math.max(blkInA.getNumRows(), blkInB.getNumRows()), Math.max(blkInA.getNumColumns(), blkInB.getNumColumns()), false);
            MatrixBlock blkOutCorr = new MatrixBlock(Math.max(blkInA.getNumRows(), blkInB.getNumRows()), Math.max(blkInA.getNumColumns(), blkInB.getNumColumns()) + 1, false);

            for (int j=0; j<smallerC; ++j) {
                int nA = (int)blkInCorrA.getValue(blkInCorrB.getNumRows() - 1, j);
                int nB = (int)blkInCorrB.getValue(blkInCorrB.getNumRows() - 1, j);

                MatrixBlock shorterMatrix = blkInA;
                MatrixBlock shorterMatrixCorr = blkInCorrA;
                int smallerR = blkInA.getNumRows();
                int largerR = blkInB.getNumRows();
                if (nB < nA) {
                    shorterMatrix = blkInB;
                    shorterMatrixCorr = blkInCorrB;
                    smallerR = blkInB.getNumRows();
                    largerR = blkInA.getNumRows();
                }

                for (int i=0; i<smallerR; ++i) {
                    int k = (int)blkInA.getValue(i, j);
                    int v = (int)blkInCorrA.getValue(i, j);
                    hashBuckets.put(k, v);

                    k = (int)blkInB.getValue(i, j);
                    v = (int)blkInCorrB.getValue(i, j);
                    hashBuckets.put(k, v);
                }

                // Todo fix bug: need taller matrix here instead of shorter
                for (int i=smallerR; i<largerR; ++i) {
                    int k = (int)shorterMatrix.getValue(i, j);
                    int v = (int)shorterMatrixCorr.getValue(i, j);
                    hashBuckets.put(k, v);
                }

                int i=0;
                for (Pair<Integer, Integer> kvPair : hashBuckets) {
                    int k = kvPair.getKey();
                    int v = kvPair.getValue();
                    blkOut.setValue(i, j, k);
                    blkOutCorr.setValue(i, j, v);
                }

                blkOutCorr.setValue(blkOutCorr.getNumRows() - 1, j, hashBuckets.size());

                hashBuckets.clear();
            }

            for (int j=smallerC; j<largerC; ++j) {
                int nHashes = (int)widerMatrix.getValue(widerMatrix.getNumRows() - 1, j);
                for (int i=0; i<nHashes; ++i) {
                    int k = (int)widerMatrix.getValue(i, j);
                    int v = (int)widerMatrixCorr.getValue(i, j);

                    blkOut.setValue(i, j, k);
                    blkOutCorr.setValue(i, j, v);
                }

                blkOutCorr.setValue(blkOutCorr.getNumRows() - 1, j, nHashes);
            }

            return new CorrMatrixBlock(blkOut, blkOutCorr);

        } else {
            throw new IllegalArgumentException("Invalid direction");
        }
    }

    @Override
    public CorrMatrixBlock intersection(CorrMatrixBlock arg0, CorrMatrixBlock arg1) {
        throw new NotImplementedException(String.format("%s intersection has not been implemented yet",
                HLLSketch.class.getSimpleName()));
    }
}
