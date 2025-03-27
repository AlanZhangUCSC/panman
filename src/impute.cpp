
#include "panmanUtils.hpp"

void panmanUtils::Tree::imputeNs(int allowedIndelDistance) {
    std::vector< std::pair< std::string, panmanUtils::NucMut > > substitutions;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > > insertions;
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > > originalNucs;
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > > wasBlockInv;

    // Make pre-order pass over the tree, building lookup tables
    fillImputationLookupTables(substitutions, insertions, originalNucs, wasBlockInv);

    // Impute all substitutions (100% success rate)
    int totalSubNs = 0;
    for (const auto& toImpute: substitutions) {
        totalSubNs += imputeSubstitution(allNodes[toImpute.first]->nucMutation, toImpute.second);
    }
    std::cout << "Imputed " << totalSubNs << "/" << totalSubNs << " SNPs/MNPs to N" << std::endl;
    
    // Attempt to impute insertions

    // Store {node ID to move : {node to move under, new mutations}} for all imputation attempts
    std::unordered_map< std::string, std::pair< panmanUtils::Node*, panmanUtils::MutationList > > toMove;
    // Must filter "insertions" to just those with Ns
    std::vector<panmanUtils::IndelPosition> insertionsWithNs;
    int insertionImputationAttempts = 0;

    // Find possible places to move nodes to for insertion imputation
    for (const auto& toImpute: insertions) {
        insertionsWithNs.clear();
        for (const auto& curInsertion: toImpute.second) {
            if (curInsertion.second > 0) {
                insertionsWithNs.push_back(curInsertion.first);
            }
        }

        // Only attempt an imputation if necessary
        if (!insertionsWithNs.empty()) {
            insertionImputationAttempts++;
            toMove[toImpute.first] = findInsertionImputationMove(
                allNodes[toImpute.first], insertionsWithNs,
                allowedIndelDistance, insertions, originalNucs, wasBlockInv);
        }
    }

    // Make all moves
    std::vector<panmanUtils::Node*> oldParents;
    std::unordered_set<panmanUtils::Node*> moved;
    for (const auto& curMove: toMove) {
        Node* curNode = allNodes[curMove.first];
        Node* newParent = curMove.second.first;
        Node* curParent = curNode->parent;
        if (newParent != nullptr) {
            // If a node was moved, any mutations calculated relative to it are no longer valid
            if (moved.find(newParent) == moved.end() && !newParent->isDescendant(moved)) {
                if (moveNode(curNode, newParent, curMove.second.second)) {
                    // This move succeeded
                    oldParents.push_back(curParent);
                    moved.emplace(curNode);
                }
            }
        }
    }

    // Compress parents with single children left over from moves
    for (const auto& curParent: oldParents) {
        if (curParent->children.size() == 1) {
            mergeNodes(curParent, curParent->children[0]);
        }
    }

    std::cout << "Moved " << moved.size() << "/" << insertionImputationAttempts << " nodes with insertions to N" << std::endl;

    // Fix depth/level attributes, post-move
    size_t numLeaves;
    size_t totalLeafDepth;
    fixLevels(root, numLeaves, totalLeafDepth);
    m_meanDepth = totalLeafDepth / numLeaves;
}

const void panmanUtils::Tree::fillImputationLookupTables( 
    std::vector< std::pair < std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {
    
    // Prepare current-state trackers
    sequence_t curSequence;
    blockExists_t blockExists;
    blockStrand_t blockStrand;
    getSequenceFromReference(curSequence, blockExists, blockStrand, root->identifier);

    // Never used, but needed so that the key is in the map
    insertions[root->identifier] = std::unordered_map< panmanUtils::IndelPosition, int32_t >();
    originalNucs[root->identifier] = std::unordered_map< panmanUtils::Coordinate, int8_t >();
    wasBlockInv[root->identifier] = std::unordered_map< uint64_t, bool >();

    for (const auto& child: root->children) {
        fillImputationLookupTablesHelper(child, substitutions, insertions, originalNucs, 
                                         wasBlockInv, curSequence, blockStrand);
    }
}

const void panmanUtils::Tree::fillImputationLookupTablesHelper(panmanUtils::Node* node, 
    std::vector< std::pair < std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv,
    sequence_t& curSequence, blockExists_t& blockStrand) {

    if (node == nullptr) return;

    fillNucleotideLookupTables(node, curSequence, substitutions, insertions, originalNucs);
    fillBlockLookupTables(node, blockStrand, wasBlockInv);

    for(auto child: node->children) {
        fillImputationLookupTablesHelper(child, substitutions, insertions, originalNucs,
                                         wasBlockInv, curSequence, blockStrand);
    }

    // Undo mutations before passing back up the tree
    for (const auto& curMut: node->nucMutation) {
        for(int i = 0; i < curMut.length(); i++) {
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);
            char originalNuc = panmanUtils::getNucleotideFromCode(originalNucs[node->identifier][curPos]);
            curPos.setSequenceBase(curSequence, originalNuc);
        }
    }

    for (const auto& curMut: node->blockMutation) {
        if (curMut.inversion) { 
            if(curMut.secondaryBlockId != -1) {
                blockStrand[curMut.primaryBlockId].second[curMut.secondaryBlockId] = 
                    !blockStrand[curMut.primaryBlockId].second[curMut.secondaryBlockId];
            } else {
                blockStrand[curMut.primaryBlockId].first = !blockStrand[curMut.primaryBlockId].first;
            }
        }
    }
}

const void panmanUtils::Tree::fillNucleotideLookupTables(panmanUtils::Node* node, sequence_t& curSequence,
    std::vector< std::pair< std::string, panmanUtils::NucMut > >& substitutions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs) {
    
    std::vector< std::pair< panmanUtils::IndelPosition, int32_t > > curNodeInsertions;
    // Will store the parent's nucleotide at all positions with an insertion, to allow reversability
    originalNucs[node->identifier] = std::unordered_map< panmanUtils::Coordinate, int8_t >();

    for (const auto& curMut: node->nucMutation) {
        int numNs = 0;
        // Handle each base in the current mutation
        for(int i = 0; i < curMut.length(); i++) {
            int8_t curNucCode = curMut.getNucCode(i);
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curMut, i);

            numNs += (curNucCode == panmanUtils::NucCode::N);

            // Mutate this position, but store the original value
            originalNucs[node->identifier][curPos] = panmanUtils::getCodeFromNucleotide(curPos.getSequenceBase(curSequence));
            curPos.setSequenceBase(curSequence, panmanUtils::getNucleotideFromCode(curNucCode));
        }

        // Save mutation if relevant
        if (curMut.isSubstitution()) {
            if (numNs > 0) {
                substitutions.emplace_back(node->identifier, curMut);
            }
        } else if (curMut.isInsertion()) {
            if (!curNodeInsertions.empty() && curNodeInsertions.back().first.mergeIndels(curMut)) {
                // Current insertion was merged with the previous one.
                curNodeInsertions.back().second += numNs;
            } else {
                // Start new insertion
                curNodeInsertions.emplace_back(panmanUtils::IndelPosition(curMut), numNs);
            }
        }
    }

    // Store curNodeInsertions as this node's insertions
    std::copy(curNodeInsertions.begin(), curNodeInsertions.end(), 
              std::inserter(insertions[node->identifier], insertions[node->identifier].begin()));
}

const void panmanUtils::Tree::fillBlockLookupTables(panmanUtils::Node* node, blockExists_t& blockStrand,
    std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {
    
    wasBlockInv[node->identifier] = std::unordered_map< uint64_t, bool >();
    for (const auto& curMut: node->blockMutation) {
        // Store original state for all deletions
        if (curMut.isDeletion()) {
            bool originalState;
            if(curMut.secondaryBlockId != -1) {
                originalState = !blockStrand[curMut.primaryBlockId].second[curMut.secondaryBlockId];
            } else {
                originalState= !blockStrand[curMut.primaryBlockId].first;
            }
            uint64_t curID = curMut.singleBlockID();
            wasBlockInv[node->identifier][curMut.singleBlockID()] = originalState;
        }
    }
}

const int panmanUtils::Tree::imputeSubstitution(std::vector<panmanUtils::NucMut>& nucMutation, const NucMut& mutToN) {
    // Get rid of the old mutation in the node's list
    std::vector<NucMut>::iterator oldIndex = std::find(nucMutation.begin(), nucMutation.end(), mutToN);
    oldIndex = nucMutation.erase(oldIndex);
    int subNs = mutToN.length();

    // Possible MNP
    if (mutToN.type() == panmanUtils::NucMutationType::NS) {
        std::vector<NucMut> snps;
        // Add non-N mutations back in (for MNPs which are partially N)
        for(int i = 0; i < mutToN.length(); i++) {
            if (mutToN.getNucCode(i) != panmanUtils::NucCode::N) {
                snps.push_back(NucMut(mutToN, i));
            }
        }
        nucMutation.insert(oldIndex, snps.begin(), snps.end());
        // These SNPs were not erased
        subNs -= snps.size();
    }
    return subNs;
}

const void panmanUtils::Tree::imputeAllSubstitutionsWithNs(std::vector<panmanUtils::NucMut>& nucMutation) {
    // Loop over vector backwards as elements will be erased
    for (int i = nucMutation.size() - 1; i >= 0; i--) {
        panmanUtils::NucMut curMut = nucMutation[i];

        if (curMut.isSubstitution()) {
            for (int i = 0; i < curMut.length(); i++) {
                if (curMut.getNucCode(i) == panmanUtils::NucCode::N) {
                    imputeSubstitution(nucMutation, curMut);
                    break;
                }
            }
        }
    }
}

const std::pair< panmanUtils::Node*, panmanUtils::MutationList > panmanUtils::Tree::findInsertionImputationMove(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    const std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {
    // Certain cases are simply impossible
    if (node == nullptr) {
        return std::make_pair(nullptr, panmanUtils::MutationList());
    }

    // Tracking best new position so far
    int bestNucImprovement = -1;
    int bestBlockImprovement = 0;
    Node* bestNewParent = nullptr;
    panmanUtils::MutationList bestNewMuts;

    std::vector<panmanUtils::IndelPosition> toImpute;
    for (const auto& curInsertion: insertions.at(node->identifier)) {
        if (curInsertion.second > 0) toImpute.push_back(curInsertion.first);
    }
    
    for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance, node,
                                                  insertions, originalNucs, wasBlockInv)) {
        panmanUtils::MutationList curNewMuts = nearby.second.concat(MutationList(node));
        curNewMuts.nucMutation = consolidateNucMutations(curNewMuts.nucMutation);
        imputeAllSubstitutionsWithNs(curNewMuts.nucMutation);
        curNewMuts.blockMutation = consolidateBlockMutations(curNewMuts.blockMutation);

        // Parsimony improvement score is the decrease in mutation count
        int nucImprovement = 0;
        int blockImprovement = node->blockMutation.size() - curNewMuts.blockMutation.size();
        for (const auto& curMut: node->nucMutation) nucImprovement += curMut.length();
        for (const auto& curMut: curNewMuts.nucMutation) nucImprovement -= curMut.length();

        if (nucImprovement > bestNucImprovement & blockImprovement >= bestBlockImprovement) {
            bestNucImprovement = nucImprovement;
            bestBlockImprovement = blockImprovement;
            bestNewParent = nearby.first;
            bestNewMuts = curNewMuts;
        }
    }

    return std::make_pair(bestNewParent, bestNewMuts);
}

const std::vector<std::pair< panmanUtils::Node*, panmanUtils::MutationList >> panmanUtils::Tree::findNearbyInsertions(
    panmanUtils::Node* node, const std::vector<panmanUtils::IndelPosition>& mutsToN, int allowedDistance, panmanUtils::Node* ignore,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::IndelPosition, int32_t > >& insertions,
    const std::unordered_map< std::string, std::unordered_map< panmanUtils::Coordinate, int8_t > >& originalNucs,
    const std::unordered_map< std::string, std::unordered_map< uint64_t, bool > >& wasBlockInv) {

    std::vector<std::pair< panmanUtils::Node*, panmanUtils::MutationList >> nearbyInsertions;

    // Bases cases: nonexistant node or node too far away
    if (node == nullptr || allowedDistance < 0) return nearbyInsertions;

    std::string curID = node->identifier;
    for (const auto& curMut: mutsToN) {
        if (insertions.at(curID).find(curMut) != insertions.at(curID).end()) {
            // Only use if this insertion has non-N nucleotides to contribute
            if (insertions.at(curID).at(curMut) < curMut.length) {
                nearbyInsertions.emplace_back(node, MutationList());
            }
            break;
        }
    }

    // Try children
    for (const auto& child: node->children) {
        if (child != ignore) {
            auto childPossibilities = findNearbyInsertions(
                child, mutsToN, allowedDistance - child->branchLength, 
                node, insertions, originalNucs, wasBlockInv);
            
            // Only bother with getting/inverting mutations if necessary
            if (!childPossibilities.empty()) {
                // Add mutations to get to child (which must be reversed)
                panmanUtils::MutationList toAdd = MutationList(child);
                toAdd.invertMutations(originalNucs.at(child->identifier), wasBlockInv.at(child->identifier));
                for (const auto& nearby: childPossibilities) {
                    nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(toAdd));
                }
            }
        }
    }
    // Try parent
    if (node->parent != ignore) {
        for (const auto& nearby: findNearbyInsertions(node->parent, mutsToN, allowedDistance - node->branchLength,
                                                      node, insertions, originalNucs, wasBlockInv)) {
            // Add mutations to get to parent
            nearbyInsertions.emplace_back(nearby.first, nearby.second.concat(MutationList(node)));
        }
    }
    return nearbyInsertions;
}

bool panmanUtils::Tree::moveNode(panmanUtils::Node* toMove, panmanUtils::Node* newParent, panmanUtils::MutationList newMuts) {
    // Avoid looping
    if (newParent->isDescendant({toMove})) return false;

    // Make dummy parent from grandparent -> dummy -> newParent
    panmanUtils::Node* dummyParent = new Node(newParent, newInternalNodeId());
    allNodes[dummyParent->identifier] = dummyParent;

    newParent->changeParent(dummyParent);
    toMove->changeParent(dummyParent);

    // newParent now has a 0-length branch from the dummy
    newParent->nucMutation.clear();
    newParent->branchLength = 0;

    // TODO: figure out how branch length works
    toMove->branchLength = 1;
    toMove->nucMutation = newMuts.nucMutation;
    toMove->blockMutation = newMuts.blockMutation;
}


void missingDataSummaryHelper(
  panmanUtils::Node* node, uint32_t& dfsIndex, sequence_t& curSequence, blockExists_t& blockExists, blockStrand_t& blockStrand,
  std::vector<std::unordered_set<panmanUtils::Coordinate>>& isMissingDataByBlock, 
  std::unordered_map<panmanUtils::Coordinate, std::vector<panmanUtils::Node*>>& missingDataByCoordinate,
  std::vector<std::pair<panmanUtils::Node*, panmanUtils::Coordinate>>& subsitutionsToMissingData
) {
  if (node == nullptr) return;

  // std::vector<std::pair<coordinate object, original nucleotide code>>
  std::vector<std::pair<panmanUtils::Coordinate, int8_t>> nucBacktrack;

  // std::vector<std::pair<coordinate object, old is missing data>>
  std::vector<std::pair<std::pair<panmanUtils::Coordinate, bool>, int32_t>> isMissingDataBacktrack;

  for (const auto& curNucMut: node->nucMutation) {
    for (int i = 0; i < curNucMut.length(); i++) {
      panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curNucMut, i);
      const auto& curPosPrimaryBlockId = curPos.primaryBlockId;
      const auto& curPosSecondaryBlockId = curPos.secondaryBlockId;
      char originalNuc = curPos.getSequenceBase(curSequence);
      int8_t originalNucCode = panmanUtils::getCodeFromNucleotide(originalNuc);
      int8_t curNucCode = curNucMut.getNucCode(i);

      if (originalNucCode == curNucCode) {
        continue;
      }

      // std::cout << curPos.primaryBlockId << ":" << curPos.secondaryBlockId << ":" << curPos.nucPosition << ":" << curPos.nucGapPosition << " " << originalNucCode << "/" << originalNuc << " -> " << curNucCode << "/" << panmanUtils::getNucleotideFromCode(curNucMut.getNucCode(i)) << std::endl;

      if (curNucCode == panmanUtils::NucCode::N) {
        subsitutionsToMissingData.emplace_back(node, curPos);
        if (isMissingDataByBlock[curPosPrimaryBlockId].find(curPos) != isMissingDataByBlock[curPosPrimaryBlockId].end()) {
          std::cerr << "Error: should not be missing data here " << curPos.primaryBlockId << ":" << curPos.secondaryBlockId << ":" << curPos.nucPosition << ":" << curPos.nucGapPosition << std::endl;
          exit(1);
        }
        isMissingDataByBlock[curPosPrimaryBlockId].insert(curPos);
        // std::cout << "Inserting missing data at " << curPos.primaryBlockId << ":" << curPos.secondaryBlockId << ":" << curPos.nucPosition << ":" << curPos.nucGapPosition << std::endl;
        isMissingDataBacktrack.emplace_back(std::make_pair(curPos, false), curPosPrimaryBlockId);
      } else if (originalNucCode == panmanUtils::NucCode::N) {
        if (isMissingDataByBlock[curPosPrimaryBlockId].find(curPos) == isMissingDataByBlock[curPosPrimaryBlockId].end()) {
          std::cerr << "Error: should be missing data here " << curPos.primaryBlockId << ":" << curPos.secondaryBlockId << ":" << curPos.nucPosition << ":" << curPos.nucGapPosition << std::endl;
          exit(1);
        }
        isMissingDataByBlock[curPosPrimaryBlockId].erase(curPos);
        // std::cout << "Erasing missing data at " << curPos.primaryBlockId << ":" << curPos.secondaryBlockId << ":" << curPos.nucPosition << ":" << curPos.nucGapPosition << std::endl;
        isMissingDataBacktrack.emplace_back(std::make_pair(curPos, true), curPosPrimaryBlockId);
      }

      nucBacktrack.emplace_back(curPos, originalNucCode);
      curPos.setSequenceBase(curSequence, panmanUtils::getNucleotideFromCode(curNucCode));
    }
  }

  //std::vector<std::tuple<primary block id, secondary block id, old on/off, old strand reversed>>
  std::vector<std::tuple<int32_t, int32_t, bool, bool>> blockBacktrack;

  for (const auto& curBlockMut : node->blockMutation) {
    int32_t primaryBlockId = curBlockMut.primaryBlockId;
    int32_t secondaryBlockId = curBlockMut.secondaryBlockId;
    bool oldStrand;
    bool oldMut;
    if (curBlockMut.isInsertion()) {
      if (curBlockMut.secondaryBlockId != -1) {
        oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
        oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
        blockExists[primaryBlockId].second[secondaryBlockId] = true;
        blockStrand[primaryBlockId].second[secondaryBlockId] = !curBlockMut.inversion;
      } else {
        oldStrand = blockStrand[primaryBlockId].first;
        oldMut = blockExists[primaryBlockId].first;
        blockExists[primaryBlockId].first = true;
        blockStrand[primaryBlockId].first = !curBlockMut.inversion;
      }
    } else if (curBlockMut.isDeletion()) {
      if (curBlockMut.secondaryBlockId != -1) {
        oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
        oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
        blockExists[primaryBlockId].second[secondaryBlockId] = false;
        blockStrand[primaryBlockId].second[secondaryBlockId] = true;
      } else {
        oldStrand = blockStrand[primaryBlockId].first;
        oldMut = blockExists[primaryBlockId].first;
        blockExists[primaryBlockId].first = false;
        blockStrand[primaryBlockId].first = true;
      }

    } else if (curBlockMut.isSimpleInversion()) {
      if (curBlockMut.secondaryBlockId != -1) {
        oldStrand = blockStrand[primaryBlockId].second[secondaryBlockId];
        oldMut = blockExists[primaryBlockId].second[secondaryBlockId];
        blockStrand[primaryBlockId].second[secondaryBlockId] = !oldStrand;
      } else {
        oldStrand = blockStrand[primaryBlockId].first;
        oldMut = blockExists[primaryBlockId].first;
        blockStrand[primaryBlockId].first = !oldStrand;
      }
    } 
    blockBacktrack.emplace_back(primaryBlockId, secondaryBlockId, oldMut, oldStrand);
  }

  size_t curNodeMissingData = 0;
  for (int i = 0; i < blockExists.size(); i++) {
    if (blockExists[i].first) {
      for (const auto& curPos: isMissingDataByBlock[i]) {
        missingDataByCoordinate[curPos].push_back(node);
        curNodeMissingData++;
      }
    }
  }

  if (node->parent != nullptr) {
    // std::cout << "DFS index: " << dfsIndex << " " << node->identifier << " parent:" << node->parent->identifier << " missing data: " << curNodeMissingData << std::endl;
  } else {
    // std::cout << "DFS index: " << dfsIndex << " " << node->identifier << std::endl;
  }

  for (auto child: node->children) {
    dfsIndex++;
    missingDataSummaryHelper(child, dfsIndex, curSequence, blockExists, blockStrand, isMissingDataByBlock, missingDataByCoordinate, subsitutionsToMissingData);
  }

  // Undo missing data before passing back up the tree
  for (const auto& [curBacktrackInfo, primaryBlockId]: isMissingDataBacktrack) {
    const auto& [curPos, oldIsMissingData] = curBacktrackInfo;
    if (oldIsMissingData) {
      if (isMissingDataByBlock[primaryBlockId].find(curPos) != isMissingDataByBlock[primaryBlockId].end()) {
        std::cerr << "Error: should not be missing data here during backtrack" << std::endl;
        exit(1);
      }
      // std::cout << "Backtracking: inserting missing data at " << curPos.primaryBlockId << ":" << curPos.secondaryBlockId << ":" << curPos.nucPosition << ":" << curPos.nucGapPosition << std::endl;
      isMissingDataByBlock[primaryBlockId].insert(curPos);
    } else {
      if (isMissingDataByBlock[primaryBlockId].find(curPos) == isMissingDataByBlock[primaryBlockId].end()) {
        std::cerr << "Error: should be missing data here during backtrack" << std::endl;
        exit(1);
      }
      // std::cout << "Backtracking: erasing missing data at " << curPos.primaryBlockId << ":" << curPos.secondaryBlockId << ":" << curPos.nucPosition << ":" << curPos.nucGapPosition << std::endl;
      isMissingDataByBlock[primaryBlockId].erase(curPos);
    }
  }

  // Undo mutations before passing back up the tree
  for (const auto& [curPos, originalNucCode]: nucBacktrack) {
    curPos.setSequenceBase(curSequence, panmanUtils::getNucleotideFromCode(originalNucCode));
  }

  // Undo block mutations before passing back up the tree
  for (const auto& [primaryBlockId, secondaryBlockId, oldMut, oldStrand]: blockBacktrack) {
    if (secondaryBlockId != -1) {
      blockStrand[primaryBlockId].second[secondaryBlockId] = oldStrand;
      blockExists[primaryBlockId].second[secondaryBlockId] = oldMut;
    } else {
      blockStrand[primaryBlockId].first = oldStrand;
      blockExists[primaryBlockId].first = oldMut;
    }
  }
}

void panmanUtils::Tree::missingDataSummary(std::ostream& fout) {
  std::cout << "Generating summary of missing data" << std::endl;

  // object keeping track of types of missing data and counts
  // Ns, none-N IUPAC ambiguous codes
  // in chains

  // Likely an unordered_map<Coordinate, std::vector<std::pair<node pointer, dfindex>>>
  // Summarize total number of Ns, total number of non-N IUPAC ambiguous codes
  // Summarize the chains, number of chains, length of chains, etc.

  uint32_t dfsIndex = 0;
  // To keep track of coorrdinates that have missing data at the current node
  std::vector<std::unordered_set<panmanUtils::Coordinate>> isMissingDataByBlock;
  // To store all missing data on the entire tree
  std::unordered_map<panmanUtils::Coordinate, std::vector<panmanUtils::Node*>> missingDataByCoordinate;
  // To store all substitutions that are to missing data
  std::vector<std::pair<panmanUtils::Node*, panmanUtils::Coordinate>> subsitutionsToMissingData;

  // Prepare current-state trackers
  sequence_t curSequence(blocks.size() + 1);
  blockExists_t blockExists(blocks.size() + 1, {false, {}});
  blockStrand_t blockStrand(blocks.size() + 1, {true, {}});
  
  for (size_t i = 0; i < blockGaps.blockPosition.size(); i++) {
    curSequence[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
    blockExists[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], false);
    blockStrand[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i], true);
  }

  int32_t maxBlockId = 0;

  // Create consensus sequence of blocks
  for(size_t i = 0; i < blocks.size(); i++) {
    isMissingDataByBlock.push_back(std::unordered_set<panmanUtils::Coordinate>());
    int32_t primaryBlockId = ((int32_t)blocks[i].primaryBlockId);
    int32_t secondaryBlockId = ((int32_t)blocks[i].secondaryBlockId);

    maxBlockId = std::max(maxBlockId, primaryBlockId);

    for(size_t j = 0; j < blocks[i].consensusSeq.size(); j++) {
      bool endFlag = false;
      for(size_t k = 0; k < 8; k++) {
        const int nucCode = (((blocks[i].consensusSeq[j]) >> (4*(7 - k))) & 15);

        if(nucCode == 0) {
          endFlag = true;
          break;
        }

        const char nucleotide = panmanUtils::getNucleotideFromCode(nucCode);

        if(secondaryBlockId != -1) {
          curSequence[primaryBlockId].second[secondaryBlockId].push_back({nucleotide, {}});
          if (nucCode == panmanUtils::NucCode::N) {
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curSequence[primaryBlockId].second[secondaryBlockId].size() - 1, -1, primaryBlockId, secondaryBlockId);
            isMissingDataByBlock[primaryBlockId].insert(curPos);
          }
        } else {
          curSequence[primaryBlockId].first.push_back({nucleotide, {}});
          if (nucCode == panmanUtils::NucCode::N) {
            panmanUtils::Coordinate curPos = panmanUtils::Coordinate(curSequence[primaryBlockId].first.size() - 1, -1, primaryBlockId, -1);
            isMissingDataByBlock[primaryBlockId].insert(curPos);
          }
        }
      }

      if(endFlag) {
        break;
      }
    }

    // End character to incorporate for gaps at the end
    if(secondaryBlockId != -1) {
        curSequence[primaryBlockId].second[secondaryBlockId].push_back({'x', {}});
    } else {
        curSequence[primaryBlockId].first.push_back({'x', {}});
    }
  }

  curSequence.resize(maxBlockId + 1);
  blockExists.resize(maxBlockId + 1);
  blockStrand.resize(maxBlockId + 1);

  // Assigning nucleotide gaps in blocks
  for(size_t i = 0; i < gaps.size(); i++) {
    int32_t primaryBId = (gaps[i].primaryBlockId);
    int32_t secondaryBId = (gaps[i].secondaryBlockId);

    for(size_t j = 0; j < gaps[i].nucPosition.size(); j++) {
      int len = gaps[i].nucGapLength[j];
      int pos = gaps[i].nucPosition[j];
      if(secondaryBId != -1) {
        curSequence[primaryBId].second[secondaryBId][pos].second.resize(len, '-');
      } else {
        curSequence[primaryBId].first[pos].second.resize(len, '-');
      }
    }
  }

  // DFS to fill missingData
  missingDataSummaryHelper(root, dfsIndex, curSequence, blockExists, blockStrand, isMissingDataByBlock, missingDataByCoordinate, subsitutionsToMissingData);


  // Print out missing data summary
  fout << "Missing data summary:" << std::endl;
  fout << "Total number of substitutions/indels to missing data: " << subsitutionsToMissingData.size() << std::endl;
  size_t totalMissingData = 0;
  for (const auto& [curPos, nodes]: missingDataByCoordinate) {
    totalMissingData += nodes.size();
  }
  fout << "Total number of missing data: " << totalMissingData << std::endl;

  std::unordered_map<std::string, size_t> missingDataCountsByNode;
  std::unordered_map<panmanUtils::Coordinate, size_t> missingDataCountsByCoordinate;
  for (const auto& [curPos, nodes]: missingDataByCoordinate) {
    for (const auto& node: nodes) {
      missingDataCountsByNode[node->identifier]++;
    }
    missingDataCountsByCoordinate[curPos] = nodes.size();
  }

  // Convert maps to vectors
  std::vector<std::pair<std::string, size_t>> missingDataCountsByNodeVec(missingDataCountsByNode.begin(), missingDataCountsByNode.end());
  std::vector<std::pair<panmanUtils::Coordinate, size_t>> missingDataCountsByCoordinateVec(missingDataCountsByCoordinate.begin(), missingDataCountsByCoordinate.end());

  // Sort vectors by count in descending order
  std::sort(missingDataCountsByNodeVec.begin(), missingDataCountsByNodeVec.end(), [](const auto& a, const auto& b) {
    return a.second > b.second;
  });

  std::sort(missingDataCountsByCoordinateVec.begin(), missingDataCountsByCoordinateVec.end(), [](const auto& a, const auto& b) {
    return a.second > b.second;
  });

  // Print top 10 nodes with the most missing data
  fout << "Top 10 nodes with the most missing data:" << std::endl;
  for (size_t i = 0; i < std::min(size_t(10), missingDataCountsByNodeVec.size()); ++i) {
    fout << "Node " << missingDataCountsByNodeVec[i].first << " has " << missingDataCountsByNodeVec[i].second << " missing data" << std::endl;
  }


  // Print top 10 coordinates with the most missing data
  fout << "Top 10 coordinates with the most missing data:" << std::endl;
  for (size_t i = 0; i < std::min(size_t(10), missingDataCountsByCoordinateVec.size()); ++i) {
    const auto& coord = missingDataCountsByCoordinateVec[i].first;
    fout << "Coordinate " << coord.primaryBlockId << ":" << coord.secondaryBlockId << ":" << coord.nucPosition << ":" << coord.nucGapPosition << " has " << missingDataCountsByCoordinateVec[i].second << " missing data" << std::endl;
  }

  // Group all the nodes for each coordinate by connected components
  std::vector<std::pair<panmanUtils::Coordinate, std::vector<std::vector<panmanUtils::Node*>>>> missingDataByCoordinateVec;
  size_t numProcessed = 0;
  for (const auto& [curPos, nodes]: missingDataByCoordinate) {
    std::unordered_set<panmanUtils::Node*> nodesSet(nodes.begin(), nodes.end());
    std::vector<std::vector<panmanUtils::Node*>> connectedComponents;
    std::unordered_set<panmanUtils::Node*> visited;
    
    // DFS function to find connected components
    std::function<void(panmanUtils::Node*, std::vector<panmanUtils::Node*>&)> dfs = 
      [&](panmanUtils::Node* node, std::vector<panmanUtils::Node*>& component) {
        if (visited.find(node) != visited.end()) {
          return;
        }
        
        visited.insert(node);
        component.push_back(node);
        
        // Check parent if it's in the original list
        if (node->parent != nullptr && 
            nodesSet.find(node->parent) != nodesSet.end()) {
          dfs(node->parent, component);
        }
        
        // Check children if they're in the original list
        for (auto& child : node->children) {
          if (nodesSet.find(child) != nodesSet.end()) {
            dfs(child, component);
          }
        }
      };
    
    // Process each node in the list
    for (auto* node : nodes) {
      if (visited.find(node) == visited.end()) {
        std::vector<panmanUtils::Node*> component;
        dfs(node, component);
        connectedComponents.push_back(component);
      }
    }
    
    missingDataByCoordinateVec.push_back({curPos, connectedComponents});
    ++numProcessed;
    std::cout << "\rProcessed " << numProcessed << " / " << missingDataByCoordinate.size() << " coordinates" << std::flush;
  }

  std::unordered_map<size_t, size_t> componentSizeCounts;
  for (const auto& [coordinate, components] : missingDataByCoordinateVec) {
    for (const auto& component : components) {
      componentSizeCounts[component.size()]++;
    }
  }

  std::vector<std::pair<size_t, size_t>> componentSizeCountsVec(componentSizeCounts.begin(), componentSizeCounts.end());
  std::sort(componentSizeCountsVec.begin(), componentSizeCountsVec.end(), [](const auto& a, const auto& b) {
    return a.first < b.first;
  });

  fout << "Connected component size distribution:" << std::endl;
  for (const auto& [size, count] : componentSizeCountsVec) {
    fout << size << " " << count << std::endl;
  }

}
