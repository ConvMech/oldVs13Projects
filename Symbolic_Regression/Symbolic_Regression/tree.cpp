#include <string.h>
#include <iostream>

using namespace ::std;

enum NODETYPE
{
	plus, minus, mult, rdiv, //restricted division 
	numTerminal, symTerminal
};

struct TreeNode {
	union {
		char nodeString[8];         // String if node is function
		double nodeDouble;          // Double if node is numeric terminal
		int nodeIndex;              // Int if node is symbolic terminal
	} nodeValue;
	TreeNode * child;
	TreeNode * sibling;
	NODETYPE kind;        // NODETYPE is a special "type field"
	int arity;
	int depth;            // Depth of this node in the tree
};

void createPlusNode(TreeNode *);
void createMinusNode(TreeNode *);
void createDoubleNode(TreeNode *, double);
void addSibling(TreeNode *, TreeNode *);
void addChild(TreeNode *, TreeNode *);
void symbolicEval(TreeNode *);

void createPlusNode(TreeNode * n) {
	strcpy_s(n->nodeValue.nodeString, "+");
	n->child = 0;
	n->sibling = 0;
	n->kind = NODETYPE::plus;
	n->arity = 2;
	n->depth = 0;
}

void createMinusNode(TreeNode * n) {
	strcpy_s(n->nodeValue.nodeString, "-");
	n->child = 0;
	n->sibling = 0;
	n->kind = NODETYPE::minus;
	n->arity = 2;
	n->depth = 0;
}

void createDoubleNode(TreeNode * n, double d) {
	n->nodeValue.nodeDouble = d;
	n->child = 0;
	n->sibling = 0;
	n->kind = numTerminal;
	n->arity = 0;
	n->depth = 0;
}

// if this node has a sibling, ask *it* to add the new one
void addSibling(TreeNode * n, TreeNode * sib) {
	if (n->sibling)
		addSibling(n->sibling, sib);
	else
		n->sibling = sib;
	sib->depth = n->depth;
}

// if this node has a child, ask the child to add a sibling
void addChild(TreeNode * n, TreeNode * ch) {
	if (n->child)
		addSibling(n->child, ch);
	else
		n->child = ch;
	ch->depth = n->depth + 1;
}

// Here a switch construction is used, since a
// node has no "self-knowledge" of what it is,
// except through its "kind" member variable   
void symbolicEval(TreeNode * n) {
	switch (n->kind) {
	case NODETYPE::plus:
		if (n->child && n->child->sibling) {
			cout << "(";
			symbolicEval(n->child);
			cout << ") + (";
			symbolicEval(n->child->sibling);
			cout << ")";
		}
		else
			cout << "+";
		break;
	case NODETYPE::minus:
		if (n->child && n->child->sibling) {
			cout << "(";
			symbolicEval(n->child);
			cout << ") - (";
			symbolicEval(n->child->sibling);
			cout << ")";
		}
		else
			cout << "-";
		break;
	case numTerminal:
		cout << n->nodeValue.nodeDouble;
		break;
	default:
		break;
	}
}

// Main.cpp
// Test driver for the TreeNode struct and related
// functions.  See Tree.h and Tree.cpp


void main() {
	TreeNode *n1, *n2, *n3, *n4, *n5;
	n1 = new TreeNode;
	n2 = new TreeNode;
	n3 = new TreeNode;
	n4 = new TreeNode;
	n5 = new TreeNode;
	createPlusNode(n1);
	createDoubleNode(n2, 2);
	createMinusNode(n3);
	createDoubleNode(n4, 3);
	createDoubleNode(n5, -1.2);
	addChild(n1, n2);
	addChild(n1, n3);
	addChild(n3, n4);
	addSibling(n4, n5);
	symbolicEval(n1);
	cout << endl;
	symbolicEval(n3);
	system("pause");
}