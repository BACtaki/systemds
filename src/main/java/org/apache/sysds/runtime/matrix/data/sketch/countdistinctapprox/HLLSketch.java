package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.lang.NotImplementedException;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.sysds.runtime.instructions.spark.data.CorrMatrixBlock;
import org.apache.sysds.runtime.matrix.data.MatrixBlock;
import org.apache.sysds.runtime.matrix.operators.CountDistinctOperator;

public class HLLSketch extends CountDistinctApproxSketch {

    private static final Log LOGGER = LogFactory.getLog(HLLSketch.class.getName());

    public HLLSketch(CountDistinctOperator op) {
        super(op);
    }

    @Override
    public Integer getScalarValue(MatrixBlock blkIn) {
        return null;
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
