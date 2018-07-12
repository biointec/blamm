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

#include <stdexcept>

#include "phylotree.h"

using namespace std;

// ============================================================================
// PHYLOGENETIC TREE (PRIVATE)
// ============================================================================

void PhylogeneticTree::recPrintTree(ostream& os, int depth) const
{
        for (int i = 0; i < depth; i++)
                cout << "\t";
        os << "-- " << branchLength << " -- " << branchName << endl;
        for (const auto& it : subTree)
                it.recPrintTree(os, depth+1);
}

float PhylogeneticTree::getTotalBranchLength() const
{
        float result = branchLength;
        for (const auto& st : subTree)
                result += st.getTotalBranchLength();

        return result;
}

void PhylogeneticTree::recScaleBranchLength(float factor)
{
        branchLength *= factor;
        for (auto& st : subTree)
                st.recScaleBranchLength(factor);
}

void PhylogeneticTree::recSetFlag(bool value)
{
        flag = value;
        for (auto& st : subTree)
                st.recSetFlag(value);
}

bool PhylogeneticTree::recFlagPaths(const set<string>& names)
{
        if (names.find(branchName) != names.end())
                flag = true;

        for (auto& st : subTree)
                if (st.recFlagPaths(names))
                        flag = true;

        return flag;
}

set<string> PhylogeneticTree::getAllNames()
{
        set<string> names;

        if (!branchName.empty())
                names.insert(branchName);

        for (auto& st : subTree) {
                set<string> stNames = st.getAllNames();
                names.insert(stNames.begin(), stNames.end());
        }

        return names;
}

void PhylogeneticTree::recUnflagRoot2LCA()
{
        flag = false;

        int numTrue = 0;
        for (const auto& st : subTree)
                if (st.flag)
                        numTrue++;

        if (numTrue < 2) {
                for (auto& st : subTree)
                        if (st.flag)
                                st.recUnflagRoot2LCA();
        }
}

float PhylogeneticTree::recGetBLS() const
{
        float value = (flag) ? branchLength : 0;

        for (const auto& st : subTree)
                value += st.recGetBLS();

        return value;
}

// ============================================================================
// PHYLOGENETIC TREE (PUBLIC)
// ============================================================================

PhylogeneticTree::PhylogeneticTree(const string& str) : flag(false), branchLength(0.0)
{
        // remove everything beyond semicolon (including semicolon)
        string input = str.substr(0, str.find_last_of(";"));

        // split input in part between () and part after ()
        size_t start = input.find_first_of("(");
        size_t end = input.find_last_of(")");

        if ((start > end) || ((start != string::npos) && (end == string::npos)))
                throw runtime_error("String does not appear to be in Newick format");

        string children = (end != string::npos) ?
                input.substr(start+1, end-start-1) : string();
        string content = (end != string::npos) ?
                input.substr(end+1) : input;

        // split content in name and length
        size_t pos = content.find_last_of(":");
        branchName = (pos != string::npos) ?
                content.substr(0, pos) : content;
        string length = (pos != string::npos) ?
                content.substr(pos+1) : string();
        if (!length.empty())
                branchLength = stof(length);

        if (children.empty())
                return;

        vector<string> subTreeStr;
        size_t prev = 0;
        for (size_t i = 0, depth = 0; i < children.size(); i++) {
                if (children[i] == '(')
                        depth++;
                if (children[i] == ')')
                        depth--;
                if ((children[i] == ',') && (depth == 0)) {
                        subTreeStr.push_back(children.substr(prev, i - prev));
                        prev = i+1;
                }
        }
        subTreeStr.push_back(children.substr(prev));

        for (const auto& it : subTreeStr)
                subTree.push_back(PhylogeneticTree(it));
}

void PhylogeneticTree::normalizeBranchLength()
{
        float total = getTotalBranchLength();
        recScaleBranchLength(1.0f / total);
}

float PhylogeneticTree::getBLS(const set<string>& names)
{
        // flag all paths from root to nodes indicated by names
        recFlagPaths(names);

        // unflag path from root to LCA
        recUnflagRoot2LCA();

        // recusively sum all branch lengths of marked nodes
        float value = recGetBLS();

        // recursively flag all nodes to false
        recSetFlag(false);

        return value;
}

ostream& operator<< (ostream& os, const PhylogeneticTree& pt)
{
        pt.recPrintTree(os);
        return os;
}
