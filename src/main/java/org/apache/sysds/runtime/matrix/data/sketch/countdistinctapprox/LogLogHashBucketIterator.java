package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.sysds.runtime.matrix.data.Pair;

import java.util.Iterator;

// Todo Parameterize
public class LogLogHashBucketIterator implements Iterator<Pair<Integer, Integer>> {
    private LogLogHashBucket hashBuckets;

    public LogLogHashBucketIterator(LogLogHashBucket hashBuckets) {
        this.hashBuckets = hashBuckets;
    }

    @Override
    public boolean hasNext() {
        return false;
    }

    @Override
    public Pair<Integer, Integer> next() {
        return null;
    }
}
