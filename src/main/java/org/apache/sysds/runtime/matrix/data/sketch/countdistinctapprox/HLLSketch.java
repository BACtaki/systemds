package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.lang.NotImplementedException;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.sysds.runtime.instructions.spark.data.CorrMatrixBlock;
import org.apache.sysds.runtime.matrix.data.MatrixBlock;
import org.apache.sysds.runtime.matrix.operators.CountDistinctOperator;
import org.apache.sysds.utils.Hash;

public class HLLSketch extends CountDistinctApproxSketch {

    private static final Log LOGGER = LogFactory.getLog(HLLSketch.class.getName());

    // This is the number used to group the input values into "buckets".
    // Each input block will contain at most 1000 < 2^10 unique values, which means that
    // 10 bits are sufficient to count all unique numbers in any given matrix block.
    // The trailing bits are used to count, which leaves 32 - 10 = 22 bits to to used to determine
    // the bucket that the hash belongs to. Each bucket stores a 32bit int, so all buckets take a total of
    // 22 * 32 = 700 bits < 100 bytes of memory.
    // TODO Explore implications of using all 10 bits; do we even need to use all of them?
    private static final int N_BITS_HASH_BUCKET_KEY = 10;

    public HLLSketch(CountDistinctOperator op) {
        super(op);
    }

    @Override
    public Integer getScalarValue(MatrixBlock blkIn) {
        HyperLogLogHashBucket hashBuckets = new HyperLogLogHashBucket(N_BITS_HASH_BUCKET_KEY);
        for (int i=0; i<blkIn.getNumRows(); ++i) {
            for (int j=0; j<blkIn.getNumColumns(); ++j) {
                int hash = Hash.hash(blkIn.getValue(i, j), this.op.getHashType());
                LOGGER.debug("Object hash (decimal)=" + hash + " Object hash (binary)=" + Integer.toBinaryString(hash));

                hashBuckets.add(hash);
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
    public MatrixBlock getMatrixValue(CorrMatrixBlock blkIn) {
        return null;
    }

    @Override
    public CorrMatrixBlock create(MatrixBlock blkIn) {
        return null;
    }

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
