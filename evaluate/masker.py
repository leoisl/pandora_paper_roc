from intervaltree import IntervalTree


class Masker:
    def __init__(self, tree: IntervalTree = None):
        self.tree = tree

    def __eq__(self, other: "Masker") -> bool:
        return self.tree == other.tree
