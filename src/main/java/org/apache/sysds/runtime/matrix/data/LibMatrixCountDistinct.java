/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

package org.apache.sysds.runtime.matrix.data;

import org.apache.commons.lang.NotImplementedException;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.sysds.api.DMLException;
import org.apache.sysds.runtime.DMLRuntimeException;
import org.apache.sysds.runtime.compress.CompressedMatrixBlock;
import org.apache.sysds.runtime.data.DenseBlock;
import org.apache.sysds.runtime.data.SparseBlock;
import org.apache.sysds.runtime.instructions.spark.data.CorrMatrixBlock;
import org.apache.sysds.runtime.matrix.data.sketch.MatrixSketch;
import org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox.HLLSketch;
import org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox.KMVSketch;
import org.apache.sysds.runtime.matrix.operators.CountDistinctOperator;
import org.apache.sysds.runtime.matrix.operators.CountDistinctOperatorTypes;
import org.apache.sysds.utils.Hash.HashType;

import java.util.HashSet;
import java.util.Set;

/**
 * This class contains various methods for counting the number of distinct values inside a MatrixBlock
 */
public class LibMatrixCountDistinct {
	private static final Log LOG = LogFactory.getLog(LibMatrixCountDistinct.class.getName());

	/**
	 * The minimum number NonZero of cells in the input before using approximate techniques for counting number of
	 * distinct values.
	 */
	public static int minimumSize = 1024;

	private LibMatrixCountDistinct() {
		// Prevent instantiation via private constructor.
	}

	/**
	 * Public method to count the number of distinct values inside a matrix. Depending on which CountDistinctOperator
	 * selected it either gets the absolute number or a estimated value.
	 * 
	 * TODO: If the MatrixBlock type is CompressedMatrix, simply read the values from the ColGroups.
	 * 
	 * @param in the input matrix to count number distinct values in
	 * @param op the selected operator to use
	 * @return the distinct count
	 */
	public static int estimateDistinctValues(MatrixBlock in, CountDistinctOperator op) {
		int res = 0;
		if((op.getOperatorType() == CountDistinctOperatorTypes.KMV || op.getOperatorType() == CountDistinctOperatorTypes.HLL) &&
			(op.getHashType() == HashType.ExpHash || op.getHashType() == HashType.StandardJava)) {
			throw new DMLException("Invalid hashing configuration using " + op.getHashType() + " and " + op.getOperatorType());
		}
		// shortcut in the simplest case.
		if(in.getLength() == 1 || in.isEmpty())
			return 1;
		else if(in.getNonZeros() < minimumSize) {
			// Just use naive implementation if the number of nonZeros values size is small.
			res = countDistinctValuesNaive(in);
		}
		else {
			switch(op.getOperatorType()) {
				case COUNT:
					res = countDistinctValuesNaive(in);
					break;
				case KMV:
					res = new KMVSketch(op).getScalarValue(in);
					break;
				case HLL:
					res = new HLLSketch(op).getScalarValue(in);
					break;
				default:
					throw new DMLException("Invalid or not implemented Estimator Type");
			}
		}

		if(res == 0)
			throw new DMLRuntimeException("Impossible estimate of distinct values");
		return res;
	}

	/**
	 * Naive implementation of counting Distinct values.
	 * 
	 * Benefit Precise, but uses memory, on the scale of inputs number of distinct values.
	 * 
	 * @param in The input matrix to count number distinct values in
	 * @return The absolute distinct count
	 */
	private static int countDistinctValuesNaive(MatrixBlock in) {
		Set<Double> distinct = new HashSet<>();
		double[] data;
		if(in.isEmpty())
			return 1;
		else if(in instanceof CompressedMatrixBlock)
			throw new NotImplementedException();

		long nonZeros = in.getNonZeros();

		if(nonZeros != -1 && nonZeros < in.getNumColumns() * in.getNumRows()) {
			distinct.add(0d);
		}

		if(in.sparseBlock != null) {
			SparseBlock sb = in.sparseBlock;

			if(in.sparseBlock.isContiguous()) {
				data = sb.values(0);
				countDistinctValuesNaive(data, distinct);
			}
			else {
				for(int i = 0; i < in.getNumRows(); i++) {
					if(!sb.isEmpty(i)) {
						data = in.sparseBlock.values(i);
						countDistinctValuesNaive(data, distinct);
					}
				}
			}
		}
		else if(in.denseBlock != null) {
			DenseBlock db = in.denseBlock;
			for(int i = 0; i <= db.numBlocks(); i++) {
				data = db.valuesAt(i);
				countDistinctValuesNaive(data, distinct);
			}
		}

		return distinct.size();
	}

	private static Set<Double> countDistinctValuesNaive(double[] valuesPart, Set<Double> distinct) {
		for(double v : valuesPart) {
			distinct.add(v);
		}
		return distinct;
	}

	public static MatrixBlock countDistinctValuesFromSketch(CorrMatrixBlock arg0, CountDistinctOperator op) {
		MatrixSketch sketch;
		if (op.getOperatorType() == CountDistinctOperatorTypes.KMV) {
			sketch = new KMVSketch(op);
		} else if (op.getOperatorType() == CountDistinctOperatorTypes.HLL) {
			sketch = new HLLSketch(op);
		} else {
			throw new NotImplementedException("Not implemented yet");
		}

		return sketch.getMatrixValue(arg0);
	}

	public static CorrMatrixBlock createSketch(MatrixBlock blkIn, CountDistinctOperator op) {

		MatrixSketch sketch;
		if (op.getOperatorType() == CountDistinctOperatorTypes.KMV) {
			sketch = new KMVSketch(op);
		} else if (op.getOperatorType() == CountDistinctOperatorTypes.HLL) {
			sketch = new HLLSketch(op);
		} else {
			throw new NotImplementedException("Not implemented yet");
		}

		return sketch.create(blkIn);
	}

	public static CorrMatrixBlock unionSketch(CorrMatrixBlock arg0, CorrMatrixBlock arg1, CountDistinctOperator op) {
		MatrixSketch sketch;
		if (op.getOperatorType() == CountDistinctOperatorTypes.KMV) {
			sketch = new KMVSketch(op);
		} else if (op.getOperatorType() == CountDistinctOperatorTypes.HLL) {
			throw new NotImplementedException("Not implemented yet");
		} else {
			throw new NotImplementedException("Not implemented yet");
		}

		return sketch.union(arg0, arg1);
	}

}
