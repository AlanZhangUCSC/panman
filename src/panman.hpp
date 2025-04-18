#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <atomic>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/task_scheduler_init.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <json/json.h>
#include "panman.capnp.h"
#include "panman.pb.h"
#include "common.hpp"

#include <mutex>

#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <kj/std/iostream.h>

namespace panmanUtils {

// 4-bit int nucleotide codes
enum NucCode {
    A = 1,
    C = 2,
    G = 4,
    T = 8,
    R = 5,
    Y = 10,
    S = 6,
    W = 9,
    K = 12,
    M = 3,
    B = 14,
    D = 13,
    H = 11,
    V = 7,
    N = 15,
    MISSING = 0
};

enum NucMutationType {
    // Nucleotide Substutution
    NS = 0,
    // Nucleotide Deletion
    ND = 1,
    // Nucleotide Insertion
    NI = 2,
    // Single Nucleotide Substitution
    NSNPS = 3,
    // Single Nucleotide Insertion
    NSNPI = 4,
    // Single Nucleotide Deletion
    NSNPD = 5,
    // None
    NNONE = 2000
};

enum BlockMutationType {
    // Block Insertion
    BI = 1,
    // Block Deletion
    BD = 0,
    // Block Inversion
    BIn = 2,
    // None
    NONE = 1000
};

// Struct for representing Nucleotide Mutation
struct NucMut {
    int32_t nucPosition;
    int32_t nucGapPosition;
    int32_t primaryBlockId;
    int32_t secondaryBlockId;
    uint8_t mutInfo;
    uint32_t nucs;

    // Default constructor
    NucMut() {}

    // Create SNP mutation for MSA (optimized for memory)
    NucMut( const std::tuple< int, int8_t, int8_t>& mutationInfo ) {
        // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
        primaryBlockId = 0;
        secondaryBlockId = -1;
        nucPosition = std::get<0>(mutationInfo);
        nucGapPosition = -1;
        mutInfo = (int)std::get<1>(mutationInfo) + (1 << 4);
        setSingleNucCode(std::get<2>(mutationInfo));
    }
    
    // Create SNP mutation
    NucMut( const std::tuple< int, int, int, int, int, int >& mutationInfo ) {
        // primaryBlockId, secondaryBlockId, pos, gapPos, type, char
        primaryBlockId = std::get<0>(mutationInfo);
        secondaryBlockId = std::get<1>(mutationInfo);
        nucPosition = std::get<2>(mutationInfo);
        nucGapPosition = std::get<3>(mutationInfo);
        mutInfo = std::get<4>(mutationInfo) + (1 << 4);
        setSingleNucCode(std::get<5>(mutationInfo));
    }

    // Create non-SNP mutations from SNP mutations at consecutive positions for MSA
    NucMut(const std::vector< std::tuple< int, int8_t, int8_t > >& mutationArray,
           int start, int end) {
        primaryBlockId = 0;
        secondaryBlockId = -1;

        mutInfo = ((end - start) << 4);
        // type
        switch(std::get<1>(mutationArray[start])) {
        case panmanUtils::NucMutationType::NSNPS:
            mutInfo += panmanUtils::NucMutationType::NS;
            break;
        case panmanUtils::NucMutationType::NSNPI:
            mutInfo += panmanUtils::NucMutationType::NI;
            break;
        case panmanUtils::NucMutationType::NSNPD:
            mutInfo += panmanUtils::NucMutationType::ND;
            break;
        case panmanUtils::NucMutationType::NS:
            mutInfo += panmanUtils::NucMutationType::NS;
            break;
        case panmanUtils::NucMutationType::NI:
            mutInfo += panmanUtils::NucMutationType::NI;
            break;
        case panmanUtils::NucMutationType::ND:
            mutInfo += panmanUtils::NucMutationType::ND;
            break;
        }

        nucPosition = (int)std::get<0>(mutationArray[start]);
        nucGapPosition = -1;

        nucs = 0;
        for(int i = start; i < end; i++) {
            addNucCode(std::get<2>(mutationArray[i]), i - start);
        }

        // if (nucPosition == 0){
        //     std::cout << "\t Writing " << nucPosition << " " << 
        //                 (int)mutInfo << " " << 
        //                 nucs << " " <<
        //                 std::endl;
        // }
    }

    // Create non-SNP mutations from SNP mutations at consecutive positions
    NucMut(const std::vector< std::tuple< int, int, int, int, int, int > >& mutationArray,
           int start, int end) {
        primaryBlockId = std::get<0>(mutationArray[start]);
        secondaryBlockId = std::get<1>(mutationArray[start]);

        mutInfo = ((end - start) << 4);
        // type
        switch(std::get<4>(mutationArray[start])) {
        case panmanUtils::NucMutationType::NSNPS:
            mutInfo += panmanUtils::NucMutationType::NS;
            break;
        case panmanUtils::NucMutationType::NSNPI:
            mutInfo += panmanUtils::NucMutationType::NI;
            break;
        case panmanUtils::NucMutationType::NSNPD:
            mutInfo += panmanUtils::NucMutationType::ND;
            break;
        case panmanUtils::NucMutationType::NS:
            mutInfo += panmanUtils::NucMutationType::NS;
            break;
        case panmanUtils::NucMutationType::NI:
            mutInfo += panmanUtils::NucMutationType::NI;
            break;
        case panmanUtils::NucMutationType::ND:
            mutInfo += panmanUtils::NucMutationType::ND;
            break;
        }

        nucPosition = std::get<2>(mutationArray[start]);
        nucGapPosition = std::get<3>(mutationArray[start]);

        nucs = 0;
        for(int i = start; i < end; i++) {
            addNucCode(std::get<5>(mutationArray[i]), i - start);
        }
    }

    // Extract mutation from protobuf nucMut object
    NucMut(panman::NucMut::Reader mutation, int64_t blockId, bool blockGapExist) {
        nucPosition = mutation.getNucPosition();
        primaryBlockId = (blockId >> 32);
        mutInfo = (mutation.getMutInfo() & 0xFF);
        nucs = (mutation.getMutInfo() >> 8);
        nucs = ((nucs) << (24 - (mutInfo >> 4)*4));

        if(blockGapExist) {
            secondaryBlockId = (blockId & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }

        if(mutation.getNucGapExist()) {
            nucGapPosition = mutation.getNucGapPosition();
        } else {
            nucGapPosition = -1;
        }
    }

    NucMut(panmanOld::nucMut mutation, int64_t blockId, bool blockGapExist) {
        nucPosition = mutation.nucposition();
        primaryBlockId = (blockId >> 32);
        mutInfo = (mutation.mutinfo() & 0xFF);
        nucs = (mutation.mutinfo() >> 8);
        nucs = ((nucs) << (24 - (mutInfo >> 4)*4));

        if(blockGapExist) {
            secondaryBlockId = (blockId & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }

        if(mutation.nucgapexist()) {
            nucGapPosition = mutation.nucgapposition();
        } else {
            nucGapPosition = -1;
        }
    }

    // Subset a SNP from an MNP
    NucMut(const NucMut& other, int i) {
        primaryBlockId = other.primaryBlockId;
        secondaryBlockId = other.secondaryBlockId;
        mutInfo = panmanUtils::NucMutationType::NSNPS + (1 << 4);
        // Extract one nucleotide, then set it as the only one
        setSingleNucCode(other.getNucCode(i));

        // If gap=-1 then increment nucPosition, otherwise increment nucGapPosition
        if (other.nucGapPosition == -1) {
            nucPosition = other.nucPosition + i;
            nucGapPosition = other.nucGapPosition;
        } else {
            nucPosition = other.nucPosition;
            nucGapPosition = other.nucGapPosition + i;
        }
    }

    // Get # of nucleotides
    int length() const {
        return (mutInfo >> 4);
    }

    // Get mutation type
    uint32_t type() const {
        return (mutInfo & 0x7);
    }

    // Get ith nucleotide code
    int getNucCode(int i) const {
        return (nucs >> (4*(5-i))) & 0xF;
    }

    // Get first nucleotide code (only nuc for NSNPX types)
    int getFirstNucCode() const {
        return getNucCode(0);
    }
    
    // Set ith nucleotide code
    void addNucCode(int8_t newNuc, int i) {
        nucs += (newNuc << (4*(5-i)));
    }

    void changeNucCode(int8_t newNuc, int i) {
        int oldCode = getNucCode(i);
        nucs -= (oldCode << (4*(5-i)));
        addNucCode(newNuc, i);
    }

    // Set to have a single nucleotide (for NSNPX types)
    void setSingleNucCode(int8_t newNuc) {
        nucs = 0;
        addNucCode(newNuc, 0);
    }

    // Is this mutation either kind of substitution?
    bool isSubstitution() const {
        return (type() == panmanUtils::NucMutationType::NSNPS
                || type() == panmanUtils::NucMutationType::NS);
    }

    // Is this mutation either kind of deletion?
    bool isDeletion() const {
        return (type() == panmanUtils::NucMutationType::NSNPD
                || type() == panmanUtils::NucMutationType::ND);
    }
    
    // Is this mutation either kind of insertion?
    bool isInsertion() const {
        return (type() == panmanUtils::NucMutationType::NSNPI
                || type() == panmanUtils::NucMutationType::NI);
    }

    bool operator==(const NucMut& other) const {
        return nucPosition == other.nucPosition &&
               nucGapPosition == other.nucGapPosition &&
               primaryBlockId == other.primaryBlockId &&
               secondaryBlockId == other.secondaryBlockId &&
               mutInfo == other.mutInfo &&
               nucs == other.nucs;
    }
};

// Struct for representing a PanMAT coordinate
struct Coordinate {
    int32_t nucPosition;
    int32_t nucGapPosition;
    int32_t primaryBlockId;
    int32_t secondaryBlockId;

    // Default constructor
    Coordinate() {}

    // Create a Coordinate by position
    Coordinate(int nucPosition, int nucGapPosition, int primaryBlockId, int secondaryBlockId) {
        this->nucPosition = nucPosition;
        this->nucGapPosition = nucGapPosition;
        this->primaryBlockId = primaryBlockId;
        this->secondaryBlockId = secondaryBlockId;
    }

    // Create a Coordinate with an offset
    Coordinate(const NucMut& nm, int offset) {
        nucPosition = nm.nucPosition;
        nucGapPosition = nm.nucGapPosition;
        primaryBlockId = nm.primaryBlockId;
        secondaryBlockId = nm.secondaryBlockId;
        moveForward(offset);
    }

    // Create a Coordinate copying a NucMut
    Coordinate(const NucMut& nm) : Coordinate(nm, 0) {}

    // Get base corresponding to this Coordinate's position within a sequence_t
    char getSequenceBase(const sequence_t& seq) const {
        if(secondaryBlockId != -1) {
            if(nucGapPosition != -1) {
                return seq[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition];
            } else {
                return seq[primaryBlockId].second[secondaryBlockId][nucPosition].first;
            }
        } else {
            if(nucGapPosition != -1) {
                return seq[primaryBlockId].first[nucPosition].second[nucGapPosition];
            } else {
                return seq[primaryBlockId].first[nucPosition].first;
            }
        }
    }

    // Set base corresponding to this Coordinate's position within a sequence_t
    void setSequenceBase(sequence_t& seq, char newNuc) const {
        if(secondaryBlockId != -1) {
            if(nucGapPosition != -1) {
                seq[primaryBlockId].second[secondaryBlockId][nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId].second[secondaryBlockId][nucPosition].first = newNuc;
            }
        } else {
            if(nucGapPosition != -1) {
                seq[primaryBlockId].first[nucPosition].second[nucGapPosition] = newNuc;
            } else {
                seq[primaryBlockId].first[nucPosition].first = newNuc;
            }
        }
    }

    // Move "offset" steps forward
    void moveForward(int offset) {
        if (nucGapPosition == -1) {
            nucPosition += offset;
        } else {
            nucGapPosition += offset;
        }
    }

    bool operator==(const Coordinate& other) const {
        return nucPosition == other.nucPosition &&
               nucGapPosition == other.nucGapPosition &&
               primaryBlockId == other.primaryBlockId &&
               secondaryBlockId == other.secondaryBlockId;
    }
};

struct IndelPosition {
    Coordinate pos;
    int32_t length;

    // Create an IndelPosition copying a NucMut
    IndelPosition(const NucMut& nm) {
        pos = Coordinate(nm);
        length = nm.length();
    }

    // Merge "other" if it comes consecutively after this IndelPosition
    // Returns whether the merge occcurs
    bool mergeIndels(const NucMut& other) {
        // Different blocks are obviously nonconsecutive
        if (pos.primaryBlockId != other.primaryBlockId || pos.secondaryBlockId != other.secondaryBlockId) {
            return false;
        }
        
        // If gap=-1 then increment nucPosition, otherwise increment nucGapPosition
        if ((pos.nucGapPosition == -1 && other.nucPosition - pos.nucPosition == length) ||
            (pos.nucGapPosition != -1 && other.nucGapPosition - pos.nucGapPosition == length)) {
            length += other.length();
            return true;
        }
        return false;
    }

    bool operator==(const IndelPosition& other) const {
        return pos == other.pos && length == other.length;
    }
};

// Struct for representing Block Mutations
struct BlockMut {
    int32_t primaryBlockId;
    int32_t secondaryBlockId;

    // Whether mutation is an insertion or deletion - Strand inversions are marked by
    // `blockMutInfo=false`, but they are not deletions
    bool  blockMutInfo;

    // Whether the block is being inverted or not. In case of insertion, whether the inserted
    // block is inverted or not
    bool inversion;

    void loadFromProtobuf(panman::Mutation::Reader mutation) {
        primaryBlockId = (mutation.getBlockId() >> 32);
        if(mutation.getBlockGapExist()) {
            secondaryBlockId = (mutation.getBlockId() & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }
        blockMutInfo = mutation.getBlockMutInfo();
        // Whether the mutation is a block inversion or not. Inversion is marked by
        // `blockMutInfo = deletion` and `inversion = true`
        inversion = mutation.getBlockInversion();
    }

    void loadFromProtobuf(panmanOld::mutation mutation) {
        primaryBlockId = (mutation.blockid() >> 32);
        if(mutation.blockgapexist()) {
            secondaryBlockId = (mutation.blockid() & 0xFFFFFFFF);
        } else {
            secondaryBlockId = -1;
        }
        blockMutInfo = mutation.blockmutinfo();
        // Whether the mutation is a block inversion or not. Inversion is marked by
        // `blockMutInfo = deletion` and `inversion = true`
        inversion = mutation.blockinversion();
    }

    BlockMut(size_t blockId, std::pair< BlockMutationType, bool > type, int secondaryBId = -1) {
        primaryBlockId = blockId;
        secondaryBlockId = secondaryBId;
        if(type.first == BlockMutationType::BI) {
            blockMutInfo = true;
        } else {
            // blockMutInfo is also set to false in the case of inversions. If the mutation
            // isn't an inversion, `blockMutInfo = false` indicates deletion
            blockMutInfo = false;
        }

        if(type.second == 1) {
            // If type.second == 1 (inversion)
            inversion = true;
        } else {
            inversion = false;
        }
    }

    BlockMut() {}

    // Readable way to check if this is an insertion
    bool isInsertion() const { return blockMutInfo; }

    // Readable way to check if this is a deletion
    bool isDeletion() const { return !blockMutInfo && !inversion; }

    // Readable way to check if this is a non-insertion inversion
    bool isSimpleInversion() const { return !blockMutInfo && inversion; }

    // Make this mutation into an insertion, keeping the same block IDs
    void convertToInsertion(bool inversion) {
        blockMutInfo = true;
        inversion = inversion;
    }

    // Make this mutation into a deletion, keeping the same block IDs
    void convertToDeletion() {
        blockMutInfo = false;
        inversion = false;
    }

    // Flip this mutation's inversion marker
    void invert() {
        inversion = !inversion;
    }

    uint64_t singleBlockID() const { 
        return (primaryBlockId << 32) + secondaryBlockId;
    }
};

// List of default blocks in the global coordinate system of the PanMAT
struct Block {
    int32_t primaryBlockId;
    int32_t secondaryBlockId;

    std::vector< uint32_t > consensusSeq;
    std::string chromosomeName;

    Block(size_t primaryBlockId, std::string seq);
    // seq is a compressed form of the sequence where each nucleotide is stored in 4 bytes
    Block(int32_t primaryBlockId, int32_t secondaryBlockId, const std::vector< uint32_t >& seq);  

    uint64_t singleBlockID() const { 
        return (primaryBlockId << 32) + secondaryBlockId;
    }
};

// List of gaps in the global coordinate system of the PanMAT
struct GapList {
    std::vector< uint32_t > nucPosition;
    int32_t primaryBlockId;
    int32_t secondaryBlockId;
    std::vector< uint32_t > nucGapLength;

};



// @DEPRECATED. To be removed when secondary block ID is removed
struct BlockGapList {
    std::vector< uint32_t > blockPosition;
    std::vector< uint32_t > blockGapLength;
};


// PanMAT tree node
class Node {
  public:
    float branchLength = 0.0;
    size_t level;
    std::string identifier;
    Node* parent;
    std::vector< Node* > children;
    std::vector< NucMut > nucMutation;
    std::vector< BlockMut > blockMutation;
    std::vector< std::string > annotations;
    bool isComMutHead = false;
    int treeIndex = -1;

    Node(std::string id, float len);
    Node(std::string id, Node* par, float len);
    // Copy another node, except for its ID, children, & annotations
    Node(Node* other, std::string id);

    // Disown/remove a child from .children
    void removeChild(Node* child) {
        children.erase(std::find(children.begin(), children.end(), child));
    }

    // Rewire pointers so node becomes a child of newParent, returning the old parent
    void changeParent(Node* newParent) {
        Node* oldParent = parent;
        // Remove node from old parent's children list
        if (oldParent != nullptr) oldParent->removeChild(this);

        // Add node to new parent's children
        parent = newParent;
        if (newParent != nullptr) newParent->children.emplace_back(this);
    }

    bool isDescendant(const std::unordered_set<Node*>& others) {
        if (parent == nullptr) {
            return false;
        } else if (others.find(parent) != others.end()) {
            return true;
        } else {
            return parent->isDescendant(others);
        }
    }
};

struct MutationList {
    std::vector<NucMut> nucMutation;
    std::vector<BlockMut> blockMutation;

    // Copy mutations from a node
    MutationList(Node* node) {
        nucMutation = node->nucMutation;
        blockMutation = node->blockMutation;       
    }

    // Copy constructor
    MutationList(const MutationList& other) {
        nucMutation = other.nucMutation;
        blockMutation = other.blockMutation;
    }

    // Default constructor (empty MutationList)
    MutationList() {}

    // Convert mutations to their exact inverse, i.e. mutations from child to parent
    void invertMutations(const std::unordered_map< Coordinate, int8_t >& originalNucs,
        const std::unordered_map< uint64_t, bool >& wasBlockInv);

    // Concatenate another MutationList to the end of this one
    MutationList concat(const MutationList& other) const {
        MutationList newMuts = *this;
        newMuts.nucMutation.insert(newMuts.nucMutation.end(), other.nucMutation.begin(), other.nucMutation.end());
        newMuts.blockMutation.insert(newMuts.blockMutation.end(), other.blockMutation.begin(), other.blockMutation.end());
        return newMuts;
    }

};

// Data structure to represent a PangenomeMAT
class Tree {
  private:
    Node* createTreeFromNewickString(std::string newick);

    // In the proto file, nodes are stored in preorder. Once the tree has been generated in
    // memory, assign mutations from the proto file to the tree nodes using preorder
    // traversal
    void assignMutationsToNodes(Node* root, size_t& currentIndex,
                                std::vector<panman::Node::Reader>& storedNode);

    void assignMutationsToNodes(Node* root, size_t& currentIndex,
                                std::vector< panmanOld::node >& nodes);

    // Get the total number of mutations of given type
    int getTotalParsimonyParallel(NucMutationType nucMutType,
                                  BlockMutationType blockMutType = NONE);
    
    void getBlockMutationsParallel();

    // Run tree traversal to extract mutations in range
    panmanUtils::Node* extractPanMATSegmentHelper(panmanUtils::Node* node,
            const std::tuple< int, int, int, int >& start,
            const std::tuple< int, int, int, int >& end, const blockStrand_t& rootBlockStrand);

    // Tree traversal for FASTA writer
    void printFASTAHelper(panmanUtils::Node* root, sequence_t& sequence,
                          blockExists_t& blockExists, blockStrand_t& blockStrand, std::ostream& fout,
                          bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
    void printFASTAHelperNew(panmanUtils::Node* root, 
                          std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                          std::vector<bool>& blockExists, 
                          std::vector<bool>& blockStrand, 
                          std::ostream& fout,
                          bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
    
    std::string printFASTAUltraFastHelper(
                          const std::vector<bool>& blockSequence,
                          std::unordered_map<int, int>& blockLengths,
                          const std::vector<panmanUtils::Node*>& nodesFromTipToRoot,  
                          std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                          std::vector<bool>& blockExists, 
                          std::vector<bool>& blockStrand, 
                          bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);

    std::string printFASTAGSUltraFastHelper(
                          const std::vector<bool>& blockSequence,
                          std::unordered_map<int, int>& blockLengths,
                          const std::vector<panmanUtils::Node*>& nodesFromTipToRoot, 
                          std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                          std::vector<bool>& blockExists, 
                          std::vector<bool>& blockStrand, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart, const std::tuple< int, int, int, int >& panMATEnd, bool allIndex);


    std::pair<std::vector<std::string>, std::vector<int>> extractSequenceHelper(
                          const std::vector<bool>& blockSequence,
                          std::unordered_map<int, int>& blockLengths,
                          const std::vector<panmanUtils::Node*>& nodesFromTipToRoot,  
                          std::vector<std::vector<std::pair<char,std::vector<char>>>>& sequence,
                          std::vector<bool>& blockExists, 
                          std::vector<bool>& blockStrand, 
                          bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
    
    std::pair<std::vector<std::string>, std::vector<int>> extractSingleSequence(panmanUtils::Node* node, bool aligned=false, bool rootSeq=false, const std::tuple<int, int, int, int> &start = {-1,-1,-1,-1}, const std::tuple<int, int, int, int>& end={-1,-1,-1,-1}, bool allIndex = false);
    
    void printSingleNodeHelper(std::vector<panmanUtils::Node*> &nodeList, int nodeListIndex, sequence_t& sequence,
        blockExists_t& blockExists, blockStrand_t& blockStrand, std::ostream& fout, bool aligned, bool rootSeq, const std::tuple< int, int, int, int >& panMATStart={-1,-1,-1,-1}, const std::tuple< int, int, int, int >& panMATEnd={-1,-1,-1,-1});

    // Merge parent and child nodes when compressing subtree
    void mergeNodes(Node* par, Node* chi);

    // Used to combine their mutations at corresponding positions when parent and child
    // nodes are combined
    std::pair< int, int > replaceMutation(std::pair<int,int> oldMutation,
                                          std::pair<int, int> newMutation);

    // Iterate through mutations and combine mutations at the same position
    std::vector< NucMut > consolidateNucMutations(const std::vector< NucMut >& nucMutation);
    // Iterate through mutations and combine mutations for the same block
    std::vector<BlockMut> consolidateBlockMutations(const std::vector<BlockMut>& blockMutation);

    // Used to confirm that consolidateNucMutations worked correctly. Can be removed in
    // production
    bool debugSimilarity(const std::vector< NucMut > array1,
                         const std::vector< NucMut > array2);

    // Compress extracted subtree by combining parent and child nodes where parent has only
    // one child
    // void compressTreeParallel(Node* node, size_t level);
    void compressTreeParallel(Node* node, size_t level, const std::set< std::string >& nodeIdsToDefinitelyInclude);

    // Used in rerooting
    void dfsExpansion(Node* node, std::vector< Node* >& vec);
    Node* transformHelper(Node* node);
    void adjustLevels(Node* node);
    // Fix .level attributes, as well as .m_maxDepth
    void fixLevels(Node* node, size_t& numLeaves, size_t& totalLeafDepth);

    // Check if tree is a polytomy
    bool hasPolytomy(Node* node);

    // Check if one PanMAT coordinate is greater than or equal to the other. Only the strand
    // information of the first block needs to be provided because if the block IDs are
    // different, the strand information does not change the result
    bool panMATCoordinateGeq(const std::tuple< int, int, int, int >& coor1,
                             const std::tuple< int, int, int, int >& coor2, bool strand);

    // Check if one PanMAT coordinate is less than or equal to the other. Only the strand
    // information of the first block needs to be provided because if the block IDs are
    // different, the strand information does not change the result
    bool panMATCoordinateLeq(const std::tuple< int, int, int, int >& coor1,
                             const std::tuple< int, int, int, int >& coor2, bool strand);

    // Functions for N imputation
    // Fill "substitutions" and "insertions" with all mutations TO N
    // "substitutions" beomces a vector of (node ID, substitution with Ns) pairs
    // "insertions" becomes a map of {node ID : {insertion position : number of Ns}}
    // "originalNucs" becomes a map of {node ID : {coordianate : nucleotide}} for original nucleotides of subsitutions/deletions
    // "wasBlockInv" becomes a map of {node ID : {block IDs : is inversion}} for original block states of deletions
    const void fillImputationLookupTables(
        std::vector< std::pair < std::string, NucMut > >& substitutions,
        std::unordered_map< std::string, std::unordered_map< IndelPosition, int32_t > >& insertions,
        std::unordered_map< std::string, std::unordered_map< Coordinate, int8_t > >& originalNucs,
        std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv);
    // Run fillImputationLookupTables for subtree with root "node"
    // Keep track of "curNucs" map with {coordinate : nucleotide} at "node", used for "originalNucs"
    // Keep track of "isInv" map with {block IDs : is inversion} at "node", used for "wasBlockInv"
    const void fillImputationLookupTablesHelper(Node* node,
        std::vector< std::pair < std::string, NucMut > >& substitutions,
        std::unordered_map< std::string, std::unordered_map< IndelPosition, int32_t > >& insertions,
        std::unordered_map< std::string, std::unordered_map< Coordinate, int8_t > >& originalNucs,
        std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv,
        sequence_t& curSequence, blockExists_t& blockStrand);
    // Fill "substitutions", "insertions", and "originalNucs" using nucleotide mutations in "node"
    // Keep track of "curNucs" map with {coordinate : nucleotide} at "node", used for "originalNucs"
    const void fillNucleotideLookupTables(Node* node, sequence_t& curSequence,
        std::vector< std::pair < std::string, NucMut > >& substitutions,
        std::unordered_map< std::string, std::unordered_map< IndelPosition, int32_t > >& insertions,
        std::unordered_map< std::string, std::unordered_map< Coordinate, int8_t > >& originalNucs);
    // Fill "wasBlockInv" using block mutations in "node"
    // Keep track of "isInv" map with {block IDs : is inversion} at "node", used for "wasBlockInv"
    const void fillBlockLookupTables(Node* node, blockExists_t& blockStrand,
        std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv);
    // Impute a specific substitution in "nucMutation", "mutToN" which mutated TO N
    // Erase mutation for maximum parsimony. Break up partially-N MNPs if needed
    // Returns the number of Ns imputed
    const int imputeSubstitution(std::vector<NucMut>& nucMutation, const NucMut& mutToN);
    // Impute all substitutions with Ns within "nucMutation"
    const void imputeAllSubstitutionsWithNs(std::vector<NucMut>& nucMutation);
    // Tries to find a similar insertion for each in "mutsToN" within "allowedDistance" branch length from "node"
    // Uses "originalNucs" and "wasBlockInv" maps to help invert mutations when necessary
    // Searches in all directions but direct descendants
    // Calculates necessary change in mutations to move there and if moving would increase parsimony
    // Returns a pair of (new parent, new mutations) for a parsimony improvement
    const std::pair< Node*, MutationList > findInsertionImputationMove(
        Node* node, const std::vector<IndelPosition>& mutsToN, int allowedDistance,
        const std::unordered_map< std::string, std::unordered_map< IndelPosition, int32_t > >& insertions,
        const std::unordered_map< std::string, std::unordered_map< Coordinate, int8_t > >& originalNucs,
        const std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv);
    // Find insertions the size/position of "mutToN" within "allowedDistance" branch length from "node"
    // Uses "originalNucs" and "wasBlockInv" maps to help invert mutations when necessary
    // Don't search down the edge to "ignore"
    // Relies on a precomputed map of nodes to insertion positions                   
    const std::vector<std::pair< Node*, MutationList >> findNearbyInsertions(
        Node* node, const std::vector<IndelPosition>& mutsToN, int allowedDistance, Node* ignore,
        const std::unordered_map< std::string, std::unordered_map< IndelPosition, int32_t > >& insertions,
        const std::unordered_map< std::string, std::unordered_map< Coordinate, int8_t > >& originalNucs,
        const std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv);

    std::string newInternalNodeId() {
        return "node_" + std::to_string(++m_currInternalNode);
    }

    size_t m_currInternalNode{ 0 };
    size_t m_numLeaves{ 0 };
    size_t m_maxDepth{ 0 };
    float m_meanDepth{ 0 };

    std::unordered_map<std::string, std::vector< std::string > > annotationsToNodes;
  public:
    Node *root;
    std::vector< Block > blocks;
    std::vector< GapList > gaps;

    // @DEPRECATED: To be removed with secondary block ID
    BlockGapList blockGaps;

    // Specifies the circular offset required to print the original sequence
    std::unordered_map< std::string, int > circularSequences;

    // Specifies the block by which the rotation algorithm rotated the sequence
    std::unordered_map< std::string, int > rotationIndexes;

    // Specifies whether sequence is inverted or not by the rotation algorithm
    std::unordered_map< std::string, bool > sequenceInverted;

    std::unordered_map< std::string, Node* > allNodes;

    Tree(const panman::Tree::Reader& mainTree);
    Tree(const panmanOld::tree& mainTree);
    Tree(std::istream& fin, FILE_TYPE ftype = FILE_TYPE::PANMAT);
    Tree(std::ifstream& fin, std::ifstream& secondFin,
         FILE_TYPE ftype = FILE_TYPE::GFA, std::string reference = "");

    // Copy blocks from current tree into new tree which is rooted at one of the internal
    // nodes of the current tree. Used in split for PanMAN
    Tree(Node* newRoot, const std::vector< Block >& b, const std::vector< GapList >& g,
         std::unordered_map< std::string, int >& cs,
         std::unordered_map< std::string, int >& ri,
         std::unordered_map< std::string, bool >& si,
         const BlockGapList& bgl);
    

    void protoMATToTree(const panman::Tree::Reader& mainTree);
    void protoMATToTree(const panmanOld::tree& mainTree);

    // Print out missing data summary
    void missingDataSummary(std::ostream& fout);
    // Impute all Ns in the Tree (meant for external use)
    void imputeNs(int allowedIndelDistance);
    // Move "toMove" to be a child of "newParent", with mutations "newMuts"
    // Return whether the move was possible without making a loop
    bool moveNode(Node* toMove, Node* newParent, MutationList newMuts);

    // Fitch Algorithm on Nucleotide mutations
    int nucFitchForwardPass(Node* node, std::unordered_map< std::string, int >& states, int refState=-1);
    int nucFitchForwardPassOpt(Node* node, std::unordered_map< std::string, int >& states);
    // Default state is used in rerooting to a tip sequence. It is used to fix the state at
    // the root
    void nucFitchBackwardPass(Node* node, std::unordered_map< std::string, int >& states,
                              int parentState, int defaultState = (1<<28));
    void nucFitchBackwardPassOpt(Node* node, std::unordered_map< std::string, int >& states,
                                 int parentState, int defaultState = (1<<28));
    void nucFitchAssignMutations(Node* node, std::unordered_map< std::string, int >& states,
                                 std::unordered_map< std::string,
                                 std::pair< panmanUtils::NucMutationType, char > >& mutations,
                                 int parentState);
    void nucFitchAssignMutationsOpt(Node* node, std::unordered_map< std::string, int >& states,
                                    std::unordered_map< std::string,
                                    std::pair< panmanUtils::NucMutationType, char > >& mutations,
                                    int parentState);

    // Sankoff algorithm on Nucleotide Mutations
    std::vector< int > nucSankoffForwardPass(Node* node, std::unordered_map< std::string, std::vector< int > >& stateSets);
    std::vector< int > nucSankoffForwardPassOpt(Node* node, std::unordered_map< std::string, std::vector< int > >& stateSets);
    void nucSankoffBackwardPass(Node* node,
                                std::unordered_map< std::string, std::vector< int > >& stateSets,
                                std::unordered_map< std::string, int >& states, int parentPtr,
                                int defaultValue = (1<<28));
    void nucSankoffBackwardPassOpt(Node* node,
                                   std::unordered_map< std::string, std::vector< int > >& stateSets,
                                   std::unordered_map< std::string, int >& states, int parentPtr,
                                   int defaultValue = (1<<28));
    void nucSankoffAssignMutations(Node* node,
                                   std::unordered_map< std::string, int >& states, std::unordered_map< std::string,
                                   std::pair< panmanUtils::NucMutationType, char > >& mutations, int parentState);
    void nucSankoffAssignMutationsOpt(Node* node,
                                      std::unordered_map< std::string, int >& states, std::unordered_map< std::string,
                                      std::pair< panmanUtils::NucMutationType, char > >& mutations, int parentState);

    // Fitch algorithm on Block Mutations
    int blockFitchForwardPassNew(Node* node,
                                 std::unordered_map< std::string, int >& states);
    void blockFitchBackwardPassNew(Node* node,
                                   std::unordered_map< std::string, int >& states, int parentState,
                                   int defaultValue = (1<<28));
    void blockFitchAssignMutationsNew(Node* node,
                                      std::unordered_map< std::string, int >& states,
                                      std::unordered_map< std::string,
                                      std::pair< panmanUtils::BlockMutationType, bool > >& mutations, int parentState);

    // Sankoff algorithm on Block Mutations
    std::vector< int > blockSankoffForwardPass(Node* node, std::unordered_map< std::string,
            std::vector< int > >& stateSets);
    void blockSankoffBackwardPass(Node* node,
                                  std::unordered_map< std::string, std::vector< int > >& stateSets,
                                  std::unordered_map< std::string, int >& states, int parentPtr,
                                  int defaultValue = (1<<28));
    void blockSankoffAssignMutations(Node* node,
                                     std::unordered_map< std::string, int >& states, std::unordered_map< std::string,
                                     std::pair< panmanUtils::BlockMutationType, bool > >& mutations, int parentState);

    // void printSummary();
    void printBlockRanges(std::ostream &out);
    void printFASTAUltraFastSubset(int64_t numPairs=100000, std::string subsetPairsFile="", bool rootSeq = false, const std::tuple<int, int, int, int> &start={-1,-1,-1,-1}, const std::tuple<int, int, int, int> &end={-1,-1,-1,-1}, bool allIndex = false);

    void printSummary(std::ostream &out);
    void printBfs(Node* node = nullptr);
    void printFASTA(std::ostream& fout, bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start={-1,-1,-1,-1}, const std::tuple<int, int, int, int> &end={-1,-1,-1,-1}, bool allIndex = false);
    void printFASTANew(std::ostream& fout, bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start={-1,-1,-1,-1}, const std::tuple<int, int, int, int> &end={-1,-1,-1,-1}, bool allIndex = false);
    void printFASTAUltraFast(std::ostream& fout, bool aligned = false, bool rootSeq = false, const std::tuple<int, int, int, int> &start={-1,-1,-1,-1}, const std::tuple<int, int, int, int> &end={-1,-1,-1,-1}, bool allIndex = false);
    void printSingleNode(std::ostream& fout, const sequence_t& sequence,
                                         const blockExists_t& blockExists, const blockStrand_t& blockStrand,
                                         std::string nodeIdentifier, std::tuple< int, int, int, int > &panMATStart, std::tuple< int, int, int, int > &panMATEnd);
    void printFASTAParallel(std::ostream& fout, bool aligned = false);
    void printMAF(std::ostream& fout);

    void printMAFNew(std::ostream& fout);
    void generateSequencesFromMAF(std::ifstream& fin, std::ofstream& fout);
    void printVCFParallel(std::string reference, std::ostream& fout);
    void printVCFParallel(panmanUtils::Node* node, std::ostream& fout);
    void extractAminoAcidTranslations(std::ostream& fout, int64_t start, int64_t end);

    // Extract PanMAT representing a segment of the genome. The start and end coordinates
    // are with respect to the root sequence. The strands of the terminal blocks in all
    // sequences are assumed to be the same as their strands in the root sequence for the
    // purpose of splitting the terminal blocks during extraction
    void extractPanMATSegment(kj::std::StdOutputStream& fout, int64_t start, int64_t end);
    void extractPanMATIndex(std::ostream& fout, int64_t start, int64_t end, std::string nodeIdentifier, bool single=true);

    Node* subtreeExtractParallel(std::vector< std::string > nodeIds, const std::set< std::string >& nodeIdsToDefinitelyInclude = {});
    // Node* subtreeExtractParallel(std::vector< std::string > nodeIds);
    void writeToFile(kj::std::StdOutputStream& fout, Node* node = nullptr);
    std::string getNewickString(Node* node);
    std::string getStringFromReference(std::string reference, bool aligned = true,
                                       bool incorporateInversions=true);
    const void getSequenceFromReference(sequence_t& sequence, blockExists_t& blockExists,
                                        blockStrand_t& blockStrand, std::string reference, bool rotateSequence = false,
                                        int* rotIndex = nullptr);

    // For each node in the tree, print mutations with respect to the root node to the
    // output file
    void printMutations(std::ostream& fout);
    void printMutationsNew(std::ostream& fout);
    void printMutationsNew(std::ostream& fout, std::string& referenceString);
    void printMutationsNew(std::ostream& fout, std::vector<std::string>& nodes, std::string& referenceString);
    void printNodePaths(std::ostream& fout);

    void getBlockSequenceFromReference(block_t& sequence, bool& blockExists,
                                       bool& blockStrand, std::string reference, int64_t primaryBlockId,
                                       int64_t secondaryBlockId);

    // Split file provided as input.
    std::pair< Tree, Tree > splitByComplexMutations(const std::string& nodeId3);

    // get unaligned global coordinate
    int32_t getUnalignedGlobalCoordinate(int32_t primaryBlockId, int32_t secondaryBlockId,
                                         int32_t pos, int32_t gapPos, const sequence_t& sequence,
                                         const blockExists_t& blockExists, const blockStrand_t& blockStrand,
                                         int circularOffset = 0, bool * check = nullptr);
    std::tuple< int, int, int, int > globalCoordinateToBlockCoordinate(
        int64_t globalCoordinate,
        const sequence_t& sequence,
        const blockExists_t& blockExists,
        const blockStrand_t& blockStrand, int64_t circularOffset = 0);

    std::string getSequenceFromVCF(std::string sequenceId, std::ifstream& fin);
    bool verifyVCFFile(std::ifstream& fin);
    void vcfToFASTA(std::ifstream& fin, std::ofstream& fout);
    void annotate(std::ifstream& fin);
    std::vector< std::string > searchByAnnotation(std::string annotation);
    void convertToGFA(std::ostream& fout);
    void convertToGFAEfficient(std::ostream& fout);
    void printFASTAFromGFA(std::ifstream& fin, std::ofstream& fout);
    void getNodesPreorder(panmanUtils::Node* root, capnp::List<panman::Node>::Builder& nodesBuilder, size_t& nodeIndex);
    size_t getGlobalCoordinate(int primaryBlockId, int secondaryBlockId, int nucPosition,
                               int nucGapPosition);

    // Transforms tree such that given node becomes child of new root
    void transform(Node* node);
    void reroot(std::string sequenceName);

    

};

// Represents complex mutations like Horizontal Gene Transfer or Recombinations
struct ComplexMutation {
    char mutationType;
    size_t treeIndex1, treeIndex2, treeIndex3;
    std::string sequenceId1, sequenceId2, sequenceId3;

    // coordinates of start in parent 1
    int32_t primaryBlockIdStart1;
    int32_t secondaryBlockIdStart1;
    int32_t nucPositionStart1;
    int32_t nucGapPositionStart1;

    // coordinates of end in parent 1
    int32_t primaryBlockIdEnd1;
    int32_t secondaryBlockIdEnd1;
    int32_t nucPositionEnd1;
    int32_t nucGapPositionEnd1;

    // coordinates of start in parent 2
    int32_t primaryBlockIdStart2;
    int32_t secondaryBlockIdStart2;
    int32_t nucPositionStart2;
    int32_t nucGapPositionStart2;

    // coordinates of end in parent 2
    int32_t primaryBlockIdEnd2;
    int32_t secondaryBlockIdEnd2;
    int32_t nucPositionEnd2;
    int32_t nucGapPositionEnd2;

    ComplexMutation(char mutType, int tIndex1, int tIndex2, int tIndex3, std::string sId1,
                    std::string sId2, std::string sId3, std::tuple< int,int,int,int > t1,
                    std::tuple< int,int,int,int > t2, std::tuple< int,int,int,int > t3,
                    std::tuple< int,int,int,int > t4) {
        mutationType = mutType;
        treeIndex1 = tIndex1;
        treeIndex2 = tIndex2;
        treeIndex3 = tIndex3;

        sequenceId1 = sId1;
        sequenceId2 = sId2;
        sequenceId3 = sId3;

        primaryBlockIdStart1 = std::get<0>(t1);
        secondaryBlockIdStart1 = std::get<1>(t1);
        nucPositionStart1 = std::get<2>(t1);
        nucGapPositionStart1 = std::get<3>(t1);

        primaryBlockIdEnd1 = std::get<0>(t2);
        secondaryBlockIdEnd1 = std::get<1>(t2);
        nucPositionEnd1 = std::get<2>(t2);
        nucGapPositionEnd1 = std::get<3>(t2);

        primaryBlockIdStart2 = std::get<0>(t3);
        secondaryBlockIdStart2 = std::get<1>(t3);
        nucPositionStart2 = std::get<2>(t3);
        nucGapPositionStart2 = std::get<3>(t3);

        primaryBlockIdEnd2 = std::get<0>(t4);
        secondaryBlockIdEnd2 = std::get<1>(t4);
        nucPositionEnd2 = std::get<2>(t4);
        nucGapPositionEnd2 = std::get<3>(t4);
    }

    ComplexMutation(panman::ComplexMutation::Reader cm) {
        mutationType = (cm.getMutationType()? 'H': 'R');
        treeIndex1 = cm.getTreeIndex1();
        treeIndex2 = cm.getTreeIndex2();
        treeIndex3 = cm.getTreeIndex3();
        sequenceId1 = cm.getSequenceId1();
        sequenceId2 = cm.getSequenceId2();
        sequenceId3 = cm.getSequenceId3();

        primaryBlockIdStart1 = (cm.getBlockIdStart1() >> 32);
        secondaryBlockIdStart1 = (cm.getBlockGapExistEnd1()?
                                  (cm.getBlockIdStart1()&(0xFFFFFFFF)): -1);
        nucPositionStart1 = cm.getNucPositionStart1();
        nucGapPositionStart1 = (cm.getNucGapExistStart1()? (cm.getNucGapPositionStart1()) : -1);

        primaryBlockIdStart2 = (cm.getBlockIdStart2() >> 32);
        secondaryBlockIdStart2 = (cm.getNucGapExistStart2()?
                                  (cm.getBlockIdStart2()&(0xFFFFFFFF)): -1);
        nucPositionStart2 = cm.getNucPositionStart2();
        nucGapPositionStart2 = (cm.getNucGapExistStart2()? (cm.getNucGapPositionStart2()) : -1);

        primaryBlockIdEnd1 = (cm.getBlockIdEnd1() >> 32);
        secondaryBlockIdEnd1 = (cm.getBlockGapExistEnd1()? (cm.getBlockIdEnd1()&(0xFFFFFFFF)): -1);
        nucPositionEnd1 = cm.getNucPositionEnd1();
        nucGapPositionEnd1 = (cm.getNucGapExistEnd1()? (cm.getNucGapPositionEnd1()) : -1);

        primaryBlockIdEnd2 = (cm.getBlockIdEnd2() >> 32);
        secondaryBlockIdEnd2 = (cm.getBlockGapExistEnd2()? (cm.getBlockIdEnd2()&(0xFFFFFFFF)): -1);
        nucPositionEnd2 = cm.getNucPositionEnd2();
        nucGapPositionEnd2 = (cm.getNucGapExistEnd2()? (cm.getNucGapPositionEnd2()) : -1);
    }

    void toCapnProto(panman::ComplexMutation::Builder& cm) {
        cm.setMutationType(mutationType == 'H');
        cm.setTreeIndex1(treeIndex1);
        cm.setTreeIndex2(treeIndex2);
        cm.setTreeIndex3(treeIndex3);
        cm.setSequenceId1(sequenceId1);
        cm.setSequenceId2(sequenceId2);
        cm.setSequenceId3(sequenceId3);

        if(secondaryBlockIdStart1 != -1) {
            cm.setBlockGapExistStart1(true);
            cm.setBlockIdStart1(((int64_t)primaryBlockIdStart1 << 32)+secondaryBlockIdStart1);
        } else {
            cm.setBlockGapExistStart1(false);
            cm.setBlockIdStart1(((int64_t)primaryBlockIdStart1 << 32));
        }
        cm.setNucPositionStart1(nucPositionStart1);

        if(nucGapPositionStart1 != -1) {
            cm.setNucGapExistStart1(true);
            cm.setNucGapPositionStart1(nucGapPositionStart1);
        }

        if(secondaryBlockIdStart2 != -1) {
            cm.setBlockGapExistStart2(true);
            cm.setBlockIdStart2(((int64_t)primaryBlockIdStart2 << 32)+secondaryBlockIdStart2);
        } else {
            cm.setBlockGapExistStart2(false);
            cm.setBlockIdStart2(((int64_t)primaryBlockIdStart2 << 32));
        }
        cm.setNucPositionStart2(nucPositionStart2);

        if(nucGapPositionStart2 != -1) {
            cm.setNucGapExistStart2(true);
            cm.setNucGapPositionStart2(nucGapPositionStart2);
        }

        if(secondaryBlockIdEnd1 != -1) {
            cm.setBlockGapExistEnd1(true);
            cm.setBlockIdEnd1(((int64_t)primaryBlockIdEnd1 << 32)+secondaryBlockIdEnd1);
        } else {
            cm.setBlockGapExistEnd1(false);
            cm.setBlockIdEnd1(((int64_t)primaryBlockIdEnd1 << 32));
        }
        cm.setNucPositionEnd1(nucPositionEnd1);

        if(nucGapPositionEnd1 != -1) {
            cm.setNucGapExistEnd1(true);
            cm.setNucGapPositionEnd1(nucGapPositionEnd1);
        }

        if(secondaryBlockIdEnd2 != -1) {
            cm.setBlockGapExistEnd2(true);
            cm.setBlockIdEnd2(((int64_t)primaryBlockIdEnd2 << 32)+secondaryBlockIdEnd2);
        } else {
            cm.setBlockGapExistEnd2(false);
            cm.setBlockIdEnd2(((int64_t)primaryBlockIdEnd2 << 32));
        }
        cm.setNucPositionEnd2(nucPositionEnd2);

        if(nucGapPositionEnd2 != -1) {
            cm.setNucGapExistEnd2(true);
            cm.setNucGapPositionEnd2(nucGapPositionEnd2);
        }

        // return cm;
    }

    ComplexMutation(panmanOld::complexMutation cm) {
        mutationType = (cm.mutationtype()? 'H': 'R');
        treeIndex1 = cm.treeindex1();
        treeIndex2 = cm.treeindex2();
        treeIndex3 = cm.treeindex3();
        sequenceId1 = cm.sequenceid1();
        sequenceId2 = cm.sequenceid2();
        sequenceId3 = cm.sequenceid3();

        primaryBlockIdStart1 = (cm.blockidstart1() >> 32);
        secondaryBlockIdStart1 = (cm.blockgapexiststart1()?
                                  (cm.blockidstart1()&(0xFFFFFFFF)): -1);
        nucPositionStart1 = cm.nucpositionstart1();
        nucGapPositionStart1 = (cm.nucgapexiststart1()? (cm.nucgappositionstart1()) : -1);

        primaryBlockIdStart2 = (cm.blockidstart2() >> 32);
        secondaryBlockIdStart2 = (cm.blockgapexiststart2()?
                                  (cm.blockidstart2()&(0xFFFFFFFF)): -1);
        nucPositionStart2 = cm.nucpositionstart2();
        nucGapPositionStart2 = (cm.nucgapexiststart2()? (cm.nucgappositionstart2()) : -1);

        primaryBlockIdEnd1 = (cm.blockidend1() >> 32);
        secondaryBlockIdEnd1 = (cm.blockgapexistend1()? (cm.blockidend1()&(0xFFFFFFFF)): -1);
        nucPositionEnd1 = cm.nucpositionend1();
        nucGapPositionEnd1 = (cm.nucgapexistend1()? (cm.nucgappositionend1()) : -1);

        primaryBlockIdEnd2 = (cm.blockidend2() >> 32);
        secondaryBlockIdEnd2 = (cm.blockgapexistend2()? (cm.blockidend2()&(0xFFFFFFFF)): -1);
        nucPositionEnd2 = cm.nucpositionend2();
        nucGapPositionEnd2 = (cm.nucgapexistend2()? (cm.nucgappositionend2()) : -1);
    }

    panmanOld::complexMutation toProtobuf() {
        panmanOld::complexMutation cm;
        cm.set_mutationtype(mutationType == 'H');
        cm.set_treeindex1(treeIndex1);
        cm.set_treeindex2(treeIndex2);
        cm.set_treeindex3(treeIndex3);
        cm.set_sequenceid1(sequenceId1);
        cm.set_sequenceid2(sequenceId2);
        cm.set_sequenceid3(sequenceId3);

        if(secondaryBlockIdStart1 != -1) {
            cm.set_blockgapexiststart1(true);
            cm.set_blockidstart1(((int64_t)primaryBlockIdStart1 << 32)+secondaryBlockIdStart1);
        } else {
            cm.set_blockgapexiststart1(false);
            cm.set_blockidstart1(((int64_t)primaryBlockIdStart1 << 32));
        }
        cm.set_nucpositionstart1(nucPositionStart1);

        if(nucGapPositionStart1 != -1) {
            cm.set_nucgapexiststart1(true);
            cm.set_nucgappositionstart1(nucGapPositionStart1);
        }

        if(secondaryBlockIdStart2 != -1) {
            cm.set_blockgapexiststart2(true);
            cm.set_blockidstart2(((int64_t)primaryBlockIdStart2 << 32)+secondaryBlockIdStart2);
        } else {
            cm.set_blockgapexiststart2(false);
            cm.set_blockidstart2(((int64_t)primaryBlockIdStart2 << 32));
        }
        cm.set_nucpositionstart2(nucPositionStart2);

        if(nucGapPositionStart2 != -1) {
            cm.set_nucgapexiststart2(true);
            cm.set_nucgappositionstart2(nucGapPositionStart2);
        }

        if(secondaryBlockIdEnd1 != -1) {
            cm.set_blockgapexistend1(true);
            cm.set_blockidend1(((int64_t)primaryBlockIdEnd1 << 32)+secondaryBlockIdEnd1);
        } else {
            cm.set_blockgapexistend1(false);
            cm.set_blockidend1(((int64_t)primaryBlockIdEnd1 << 32));
        }
        cm.set_nucpositionend1(nucPositionEnd1);

        if(nucGapPositionEnd1 != -1) {
            cm.set_nucgapexistend1(true);
            cm.set_nucgappositionend1(nucGapPositionEnd1);
        }

        if(secondaryBlockIdEnd2 != -1) {
            cm.set_blockgapexistend2(true);
            cm.set_blockidend2(((int64_t)primaryBlockIdEnd2 << 32)+secondaryBlockIdEnd2);
        } else {
            cm.set_blockgapexistend2(false);
            cm.set_blockidend2(((int64_t)primaryBlockIdEnd2 << 32));
        }
        cm.set_nucpositionend2(nucPositionEnd2);

        if(nucGapPositionEnd2 != -1) {
            cm.set_nucgapexistend2(true);
            cm.set_nucgappositionend2(nucGapPositionEnd2);
        }

        return cm;
    }

};

// Data structure to represent PanMAN
class TreeGroup {
  public:
    // List of PanMATs in PanMAN
    std::vector< Tree > trees;
    // List of complex mutations linking PanMATs
    std::vector< ComplexMutation > complexMutations;

    TreeGroup(std::istream& fin, bool isOld = false);
    // List of PanMAT files and a file with all the complex mutations relating these files
    TreeGroup(std::vector< std::ifstream >& treeFiles, std::ifstream& mutationFile);
    TreeGroup(std::vector< Tree* >& t);
    TreeGroup(std::vector< Tree* >& tg, std::ifstream& mutationFile);

    TreeGroup* subnetworkExtract(std::unordered_map< int, std::vector< std::string > >& nodeIds);

    void printFASTA(std::ofstream& fout, bool rootSeq = false);
    void writeToFile(kj::std::StdOutputStream& fout);
    void printComplexMutations(std::ostream& fout);
};

};

namespace std {
    template <>
    struct hash<panmanUtils::Coordinate> {
        size_t operator()(const panmanUtils::Coordinate& coord) const {
            size_t h1 = std::hash<int32_t>{}(coord.nucPosition);
            size_t h2 = std::hash<int32_t>{}(coord.nucGapPosition);
            size_t h3 = std::hash<int32_t>{}(coord.primaryBlockId);
            size_t h4 = std::hash<int32_t>{}(coord.secondaryBlockId);
            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
        }
    };
    
    template <>
    struct hash<panmanUtils::IndelPosition> {
        size_t operator()(const panmanUtils::IndelPosition& indelPos) const {
            size_t h1 = std::hash<panmanUtils::Coordinate>{}(indelPos.pos);
            size_t h2 = std::hash<int32_t>{}(indelPos.length);
            return h1 ^ (h2 << 4);
        }
    };
}
