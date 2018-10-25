import numpy as np

class DisjointSetForest(object):
    """Efficient union find implementation for constructing the linkage matrix

    Parameters
    ----------

    size : int
        The total size of the set of objects to
        track via the union find structure.

    Attributes
    ----------

    is_component : array of bool; shape (size, 1)
        Array specifying whether each element of the
        set is the root node, or identifier for
        a component.
        
    Reference
    ----------
    
    https://en.wikipedia.org/wiki/Disjoint-set_data_structure
    It provides near-constant-time operations (bounded by the inverse Ackermann function) 
    to add new sets, to merge existing sets, and to determine whether elements are in the same set.
    """
    def __init__(self, N):
        self.parent = np.array(range(2*N-1))
        self.parent = self.parent.astype(int)
        self.next_label = N 
        self.size = np.hstack((np.ones(N),  # initial all are singleton
                                   np.zeros(N-1)))  # for cluster that are formed later

    def union(self, m, n):
        self.size[self.next_label] = self.size[m] + self.size[n]
        self.parent[m] = self.next_label
        self.parent[n] = self.next_label
        # self.size[self.next_label] = self.size[m] + self.size[n]
        self.next_label += 1
        return

    def fast_find(self, n):
        p = n
        while p != self.parent[p]:
            p = self.parent[p]

        c = n
        while c != p:
            old_parent = self.parent[c]
            self.parent[c] = p
            c = old_parent
        return p