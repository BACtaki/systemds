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

import java.util.BitSet;

/**
 * TODO
 */
public final class LongBitSet {
    private BitSet left;
    private BitSet right;

    public LongBitSet() {
        left = new BitSet();
        right = new BitSet();
    }

    public LongBitSet(long nbits) {
        throw new NotImplementedException("Constructor not implemented yet");
    }

    public boolean get(long bitIndex) {
        if (bitIndex > Integer.MAX_VALUE) {
            int diff = (int) bitIndex - Integer.MAX_VALUE;
            return left.get(diff);
        } else {
            return right.get((int) bitIndex);
        }
    }

    public LongBitSet get(int fromIndex, int toIndex) {
        throw new NotImplementedException("Method not implemented yet");
    }

    public void set(long bitIndex) {
        if (bitIndex > Integer.MAX_VALUE) {
            int diff = (int) bitIndex - Integer.MAX_VALUE;
            left.set(diff);
        } else {
            right.set((int) bitIndex);
        }
    }

    public void set(long bitIndex, boolean value) {
        throw new NotImplementedException("Method not implemented yet");
    }

    public void set(long fromIndex, long toIndex) {
        throw new NotImplementedException("Method not implemented yet");
    }

    public void set(long fromIndex, long toIndex, boolean value) {
        throw new NotImplementedException("Method not implemented yet");
    }
}
