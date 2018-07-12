/***************************************************************************
 *   Copyright (C) 2017-2018 Jan Fostier (jan.fostier@ugent.be)            *
 *   This file is part of Blamm                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef PHYLOTREE_H
#define PHYLOTREE_H

#include <string>
#include <vector>
#include <iostream>
#include <set>

// ============================================================================
// PHYLOGENETIC TREE
// ============================================================================

class PhylogeneticTree {

private:
        bool flag;                              // aux. flag for LCA computation
        float branchLength;                     // branch length
        std::string branchName;                 // branch name
        std::vector<PhylogeneticTree> subTree;  // subtree

        /**
         * Recursively print the phylogenetic tree
         * @param os Output stream
         * @param depth Recursive depth
         */
        void recPrintTree(std::ostream& os, int depth = 0) const;

        /**
         * Get the total branch length of the tree
         * @return The total branch length of the tree
         */
        float getTotalBranchLength() const;

        /**
         * Recursively scale the branch length with some proportionality factor
         * @param factor proportionality factor
         */
        void recScaleBranchLength(float factor);

        /**
         * Recursively set the flag of the tree to value
         * @param value Target value
         */
        void recSetFlag(bool value);

        /**
         * Given a set of node names, recursively flag all paths to those nodes
         * @param names Set of node names
         * @return true if at least one node was flagged
         */
        bool recFlagPaths(const std::set<std::string>& names);

        /**
         * Unflag the path from root to lowest common ancestor (LCA), including
         * the LCA itself as the LCA itself does not contribute to the BLS
         */
        void recUnflagRoot2LCA();

        /**
         * Recursively sum branch lengths for flagged nodes
         * @return The branch length score (BLS)
         */
        float recGetBLS() const;

public:
        /**
         * Default constructor
         */
        PhylogeneticTree() : flag(false), branchLength(0.0) {};

        /**
         * Create a phylogentic tree from Newick format
         * @param str String with tree description in Newick tree format
         */
        PhylogeneticTree(const std::string& str);

        /**
         * Normalize the branch lengths of the tree
         */
        void normalizeBranchLength();

        /**
         * Given a set of names, compute the BLS
         * @param names Set of node names (leafs or internal nodes)
         * @return The branch length score
         */
        float getBLS(const std::set<std::string>& names);

        /**
         * Get all branch names in the tree
         * @return Set of branch names
         */
        std::set<std::string> getAllNames();

        /**
         * operator<< overloading
         * @param os Output stream
         * @param pt PhylogeneticTree
         */
        friend std::ostream& operator<< (std::ostream& os, const PhylogeneticTree& pt);
};

#endif
