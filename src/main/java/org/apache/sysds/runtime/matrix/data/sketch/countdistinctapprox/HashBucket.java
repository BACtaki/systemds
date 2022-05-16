package org.apache.sysds.runtime.matrix.data.sketch.countdistinctapprox;

import org.apache.commons.lang.NotImplementedException;

import java.util.TreeSet;

public class HashBucket {
    private TreeSet<Integer> sortedHashes;

    public void addToBucket() {
        throw new NotImplementedException("");
    }


    public int getSmallestKPercent(int k) {
        throw new NotImplementedException("");
    }
}
