package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.sysds.runtime.matrix.data.Pair;

import java.util.HashMap;
import java.util.Iterator;

// Todo Parameterize
class HashMapPairIterator<T> implements Iterator<Pair<T, T>> {
    // Todo logger

    private HashMap<T, T> hashBuckets;
    private Iterator<T> hashBucketsKeysIterator;

    public HashMapPairIterator(HashMap<T, T> hashBuckets) {
        this.hashBuckets = hashBuckets;
        this.hashBucketsKeysIterator = this.hashBuckets.keySet().iterator();
    }

    @Override
    public boolean hasNext() {
        return this.hashBucketsKeysIterator.hasNext();
    }

    @Override
    public Pair<T, T> next() {
        if (!this.hasNext()) {
            throw new IndexOutOfBoundsException(this.getClass().getSimpleName() + " cannot iterate over empty set");
        }

        T key = hashBucketsKeysIterator.next();
        T value = this.hashBuckets.get(key);

        return new Pair<>(key, value);
    }
}
