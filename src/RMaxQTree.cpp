/*
 * RMaxQTree.cpp
 *
 * Created Oct 7th 2017
 *	Author: Anna Kuosmanen
 *
 */

#include "RMaxQTree.h"

#include <vector>
#include <limits>
#include <math.h>
#include <iostream>

//i_type negative_infinity = - std::numeric_limits<i_type>::infinity();
//i_type negative_infinity = -10000000; // TODO understand this
const i_type negative_infinity = std::numeric_limits<i_type>::min();

void RMaxQTree::init(i_type node, i_type b, i_type e, i_type *keys) {
	// leaf
	if (b == e)
		tree[node].setValues(keys[b],-1, negative_infinity);
	else {
		// split
		init(2 * node, b, (b + e) / 2, keys);
		init(2 * node + 1, (b + e) / 2 + 1, e, keys);
		// propagate up
		if (tree[2 * node].Cj >= tree[2 * node + 1].Cj)
			tree[node] = tree[2 * node];
		else
			tree[node] = tree[2 * node + 1];
	}
}
	
	
// Called by public function update with only two parameters
void RMaxQTree::updateTree(i_type key, i_type j, i_type Cj, i_type node, i_type b, i_type e) {
	// Ended in a leaf
	if (b == e) {
		if(tree[node].key == key) {
			tree[node].j = j;
			tree[node].Cj = Cj;
		}
	}
	// Search
	else {	
		i_type mid = (b + e) / 2;
		if (key <= keys[mid])
			updateTree(key, j, Cj, 2 * node, b, mid);
		else
			updateTree(key, j, Cj, 2 * node + 1, mid + 1, e);
		// And propagate back up
		if (tree[2 * node].Cj >= tree[2 * node + 1].Cj)
			tree[node] = tree[2 * node];
		else
			tree[node] = tree[2 * node + 1];
	}
}

// Called by public function query with only two parameters
// Returns pair (j,C[j])
// Note here that i and j are values of _keys_, whereas b and e are indexes in the keys array!
std::pair<i_type,i_type> RMaxQTree::queryTree(i_type i, i_type j, i_type node, i_type b, i_type e) {
	// bad i_typeerval
	if (i > keys[e] || j < keys[b])
		return std::make_pair(-1,negative_infinity);
	// good i_typeerval
	if (keys[b] >= i && keys[e] <= j)
		return std::make_pair(tree[node].j,tree[node].Cj);
	// check left and right subtree
	std::pair<i_type,i_type> left = queryTree(i, j, 2 * node, b, (b + e) / 2);
	std::pair<i_type,i_type> right = queryTree(i, j, 2 * node + 1, (b + e) / 2 + 1, e);
	if (left.second == negative_infinity)
		return right;
	if (right.second == negative_infinity)
		return left;
	if (left.second >= right.second)
		return left;
	return right;
}
	

// Empty constructor for creating arrays
RMaxQTree::RMaxQTree() {}

RMaxQTree::~RMaxQTree() {
	delete [] this->tree;
}

// For filling the empty RMaxQTrees
void RMaxQTree::fillRMaxQTree(i_type *keys, i_type keyLen) {
	this->keyLen = keyLen;
	this->keys = keys;
	this->treeLen = 2 << (i_type)ceil(log2(keyLen));
	this->tree = new TreeNode[treeLen];
	init(1, 0, keyLen - 1, keys);	
}	

RMaxQTree::RMaxQTree(i_type *keys, i_type keyLen) {
	this->keyLen = keyLen;
	this->keys = keys;
	this->treeLen = 2 << (i_type)ceil(log2(keyLen));
	this->tree = new TreeNode[treeLen];
	init(1, 0, keyLen - 1, keys);
}

void RMaxQTree::update(i_type key, i_type j, i_type Cj) {
	// Start with node 1, go over the whole tree till find the key
	this->updateTree(key, j, Cj, 1, 0, keyLen - 1);
}
	
std::pair<i_type,i_type> RMaxQTree::query(i_type start, i_type end) {
	// Go over the whole tree till find the query
	return this->queryTree(start, end, 1, 0, keyLen - 1);
}
 
