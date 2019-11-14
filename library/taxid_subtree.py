"""
Randomly permutes a fasta file and a taxid file.  Can split the data into
training, validation, and testing partitions.

:Authors:
        Jacob Porter <jsporter@vt.edu>
"""


def get_subtree(tree, root, taxid_list, childless):
    """
    Perform a depth first search of the taxonomic tree
    starting at root.  Add taxids to the taxid_list.
    """
    children = tree[root]
    if children or childless:
        taxid_list.append(root)
    for child in children:
        get_subtree(tree, child, taxid_list, childless)
