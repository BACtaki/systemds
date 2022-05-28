package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.lang.NotImplementedException;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.sysds.common.Types;
import org.apache.sysds.runtime.instructions.spark.data.CorrMatrixBlock;
import org.apache.sysds.runtime.matrix.data.MatrixBlock;
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
        int M = blkIn.getNumRows(), N = blkIn.getNumColumns();
        if (this.op.getDirection() == Types.Direction.RowCol) {
            // Todo overwrite input instead of allocating new memory
//            MatrixBlock blkOut = new MatrixBlock(blkIn);

            // M x N -> M x 1
//            MatrixBlock blkOut = new MatrixBlock(blkIn.getNumRows(), 2, false);
            MatrixBlock blkOut = new MatrixBlock((int)Math.pow(2, K), 2, false);
            MatrixBlock blkOutCorr = new MatrixBlock(1, 1, false);

            // Assume blkOut is wide enough to hold hash_bucket (key, value)
            // pairs, i.e. blkOut comprises at least 2 rows/columns

            for (int i=0; i<M; ++i) {
                for (int j=0; j<N; ++j) {
                    int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                    LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));
                    hashBuckets.addHash(hash);
                }
            }

            for ()
            hashBuckets.serialize(blkOut, op.getDirection());

            blkOutCorr.setValue(0, 0, hashBuckets.size());
            return new CorrMatrixBlock(blkOut, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Row) {
            // Todo overwrite input instead of allocating new memory
//            MatrixBlock blkOut = new MatrixBlock(blkIn);

            // M x N -> M x 1
            MatrixBlock blkOut = new MatrixBlock(blkIn.getNumRows(), 2, false);
//            MatrixBlock blkOut = new MatrixBlock((int)Math.pow(2, K), 2, false);
            MatrixBlock blkOutCorr = new MatrixBlock(blkIn.getNumRows(), 1, false);

            // Assume blkOut is wide enough to hold hash_bucket (key, value)
            // pairs, i.e. blkOut comprises at least 2 rows/columns

            for (int i=0; i<M; ++i) {
                for (int j=0; j<N; ++j) {
                    int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                    LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));
                    hashBuckets.addHash(hash);
                }

                // New row
                blkOutCorr.setValue(i, 0, hashBuckets.size());
                hashBuckets.clear();
            }

            hashBuckets.serialize(blkOut, op.getDirection());
            return new CorrMatrixBlock(blkOut, blkOutCorr);

        } else if (this.op.getDirection() == Types.Direction.Col) {
            MatrixBlock blkOut = new MatrixBlock(2, blkIn.getNumColumns(), false);
//            MatrixBlock blkOut = new MatrixBlock(2, (int)Math.pow(2, K), false);
            MatrixBlock blkOutCorr = new MatrixBlock(1, blkIn.getNumColumns(), false);

            for (int j=0; j<N; ++j) {
                for (int i=0; i<M; ++i) {
                    int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                    LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));
                    hashBuckets.addHash(hash);
                }

                // New column
                blkOutCorr.setValue(0, j, hashBuckets.size());
                hashBuckets.clear();
            }

            hashBuckets.serialize(blkOut, op.getDirection());
            return new CorrMatrixBlock(blkOut, blkOutCorr);

        } else {
            throw new IllegalArgumentException("Unrecognized direction");
        }
    }

//    private MatrixBlock sliceMatrixBlockByIndexDirection(MatrixBlock blkIn, int idx) {
//        MatrixBlock blkInSlice;
//        if (op.getDirection().isRow()) {
//            blkInSlice = blkIn.slice(idx, idx);
//        } else if (op.getDirection().isCol()) {
//            blkInSlice = blkIn.slice(0, blkIn.getNumRows() - 1, idx, idx);
//        } else {
//            blkInSlice = blkIn;
//        }
//
//        return blkInSlice;
//    }

    @Override
    public CorrMatrixBlock union(CorrMatrixBlock arg0, CorrMatrixBlock arg1) {
        return null;
    }

    @Override
    public CorrMatrixBlock intersection(CorrMatrixBlock arg0, CorrMatrixBlock arg1) {
        throw new NotImplementedException(String.format("%s intersection has not been implemented yet",
                HLLSketch.class.getSimpleName()));
    }
}
