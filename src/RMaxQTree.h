/*
 * RMaxQTree.h
 *
 * Created Oct 7th 2017
 *	Author: Anna Kuosmanen
 *
 * A Range Maximum Query Tree.
 *
 * A static implementation. All the keys must be given at the construction.
 */

#ifndef RMaxQTREE_H_
#define RMaxQTREE_H_

#include <utility>

typedef long long i_type;
class TreeNode {
public:
	i_type key; // key
	i_type Cj; // The best coverage
	i_type j; // The Tuple index that the best coverage corresponds to

	// Dummy needed for initializing the array
	TreeNode() {
	}
	
	TreeNode(i_type key, i_type j, i_type Cj) {
		this->key = key;
		this->Cj = Cj;
		this->j = j;
	}

	void setValues(i_type key, i_type j, i_type Cj) {
		this->key = key;
		this->Cj = Cj;
		this->j = j;
	}
};

class RMaxQTree {

private:
	TreeNode *tree;
	i_type *keys;
	i_type treeLen, keyLen;

	void init(i_type node, i_type b, i_type e, i_type *keys);

	// Recursive helper function
	void updateTree(i_type key, i_type j, i_type Cj, i_type node, i_type b, i_type e);

	// Recursive helper function
	std::pair<i_type,i_type> queryTree(i_type i, i_type j, i_type node, i_type b, i_type e);

public:

	// Empty constructor for creating arrays
	RMaxQTree();

	~RMaxQTree();

	// For filling the empty RMaxQTrees
	void fillRMaxQTree(i_type *keys, i_type keyLen);

	RMaxQTree(i_type *keys, i_type keyLen);

	// Update the node with key "key" with the given values
	void update(i_type key, i_type j, i_type Cj);

	// Return pair (j,C[j])
	std::pair<i_type,i_type> query(i_type start, i_type end);
};

#endif /* _RMaxQTREE_H_ */
