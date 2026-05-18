# object class for construction of rooted and unrooted binary trees.
#
# lawrence david - ldavid@mit.edu - 2006.

import re
import copy
import random
import math
import sys
import pdb

from numpy import *

import rand_node as node
import rand_branch as branch


class multitree:

    def __init__(self):
        self.build_queue = []
        self.a_node = None
        self.a_branch = None
        self.root = None
        self.root_branch = None
        self.num_dups = 0
        self.num_losses = 0
        self.num_xfers = 0
        self.event_queue = []
        self.node_dict = dict()
        self.mapping_hash = dict()
        self.loss_dict = {}
        self.lca_dict = {}
        self.rec_dict = {}
        self.leaf_count = 0
        self.species_count = {}
        self.internal_node_list = []
        self.leaf_node_list = []
        self.branch_list = []
        self.liks = None
        self.migration_matrix = None
        self.filters = None
        self.thresh_dict = {}
        self.base_colors = None

    def __repr__(self):
        if self.root is None:
            return "not rooted"
        else:
            newickString = ""
            return self.root.treePrint(newickString)

    def PrintLabeledBoots(self):

        all_liks = array([i.sum_lik for i in self.internal_node_list])
        self.liks = all_liks

        if self.root is None:
            return "not rooted"
        else:
            newickString = ""
            return self.root.BootPrint(newickString, 1)

    def GetSpeciesDict(self):
        species_count = {}

        for i in self.node_dict.values():
            if i.isLeaf():
                if i.species not in species_count:
                    species_count[i.species] = 1
                else:
                    species_count[i.species] += 1
        self.species_count = species_count
        return self.species_count

    def build(self, newick):

        p = re.compile(r'\)[\d\.]+:')
        newick = p.sub('):', newick)
        self.RecursiveBuild(newick)

    def RecursiveBuild(self, newick):

        colon_match1 = re.findall(':', newick)
        colon_match2 = re.findall(r':[^\)]*:[^\)]*:', newick)

        if len(colon_match1) == 1:

            root_node = self.build_queue[0]

            for kid_branches in root_node.branch_list:
                root_node.child_branches.append(kid_branches)

            self.root = root_node
            self.root.imposeHierarchy()
            self.labelSubtrees()
            self.Get_Subnodes()
            self.root.Fill_Node_Dict()

        elif len(colon_match1) == 2:

            regex = re.search(r'\([^\(\)]*\)', newick).span()
            clade = newick[regex[0] + 1:regex[1] - 1]

            children = re.split(',', clade)
            son = children[0]
            daughter = children[1]

            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)

            center_node = node1.unite(node2)

            dummy_node = self.BuildNode("dummy_node:0.1")
            dummy_node.branch_list = []
            self.leaf_count -= 1
            dummy_node.addBranch(0.01)
            dummy_node.myBranch.addNode(center_node)
            self.a_branch = dummy_node.myBranch
            self.a_node = dummy_node
            center_node.UnrootedLeaving()

        elif len(colon_match1) == 3 and len(colon_match2) > 0:

            regex = re.search(r'\([^,]*,', newick).span()
            first_leaf = newick[regex[0] + 1:regex[1] - 1]

            newick = newick.replace(first_leaf + ',', '')

            first_node = self.BuildNode(first_leaf)

            self.build_queue.append(first_node)

            regex = re.search(r'\(.*,[^\)]*\);', newick).span()
            clade = newick[regex[0] + 1:regex[1] - 2]

            children = re.split(',', clade)
            son = children[0]
            daughter = children[1]

            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)
            center_node = node1.unite(node2)

            last_node = self.build_queue[0]
            regex = re.search(r':[^,]*,', newick).span()
            last_node.myBranch.addNode(center_node)

            self.a_branch = last_node.myBranch

            if len(self.a_branch.ends[0].branch_list) == 3:
                self.a_node = self.a_branch.ends[0]
                self.a_branch.ends[0].UnrootedLeaving()
            else:
                self.a_node = self.a_branch.ends[1]
                self.a_branch.ends[1].UnrootedLeaving()

        else:

            regex = re.search(r'\([^\(\)]*\):', newick).span()
            clade = newick[regex[0] + 1:regex[1] - 2]

            children = re.split(',', clade)
            son = children[0]
            daughter = children[1]

            node1 = self.BuildNode(son)
            node2 = self.BuildNode(daughter)
            center_node = node1.unite(node2)

            new_newick = newick[0:regex[0]]
            new_newick += center_node.name
            new_newick += newick[regex[1] - 1:]

            self.build_queue.append(center_node)

            self.RecursiveBuild(new_newick)

    def BuildNode(this_tree, node_string):
        splitSon = re.split(':', node_string)

        node1match = False
        for i in range(len(this_tree.build_queue)):
            if this_tree.build_queue[i].name == splitSon[0]:
                node1 = this_tree.build_queue[i]
                this_tree.build_queue.remove(node1)
                node1match = True
                break
        if node1match is not True:
            node1 = node.node(splitSon[0], this_tree)
            this_tree.leaf_count += 1

        node1.addBranch(splitSon[1])

        return node1

    def labelSubtrees(self):
        if self.root is None:
            print("you're trying to label an unrooted tree")
        else:
            self.root.subtreeLabel()

    def relabelSubtrees(self):
        if self.root is None:
            print("you're trying to label an unrooted tree")
        else:
            self.root.subtreeReLabel()

    def Get_Subnodes(self):
        self.root.Find_Subnodes()

    def rootify(self, center_branch):

        self.root_branch = center_branch

        center_name = center_branch.ends[0].name
        center_name += "-" + center_branch.ends[1].name
        center_node = node.node(center_name, self)
        self.root = center_node

        child1 = branch.branch(center_branch.length / 2)
        child1.addNode(center_node)
        child1.addNode(center_branch.ends[0])
        child2 = branch.branch(center_branch.length / 2)
        child2.addNode(center_node)
        child2.addNode(center_branch.ends[1])
        center_node.child_branches.append(child1)
        center_node.child_branches.append(child2)

        for kids in center_branch.ends:
            kids.branch_list.remove(center_branch)

        center_node.imposeHierarchy()
        self.labelSubtrees()
        self.Get_Subnodes()
        self.root.Fill_Node_Dict()

    def LeafShuffle(this_tree):

        # get list of leaves
        leaves = this_tree.leaf_node_list

        # assign each a random number
        shuffle_dict = {}
        for leaf in leaves:
            rand_ind = random.random()
            shuffle_dict[rand_ind] = leaf

        keys = sorted(shuffle_dict.keys())

        # swap in new species assignment
        for ind in range(len(leaves)):
            donor = leaves[ind]
            accep = shuffle_dict[keys[ind]]
            accep.perturb_species = donor.species

        return
