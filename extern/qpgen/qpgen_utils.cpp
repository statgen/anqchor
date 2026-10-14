#include "qpgen_utils.h"

#include <map>

// identify overlapping IDS from list1 and list2, and store sorted overlapping IDs in overlapping_ids
int32_t identify_overlapping_ids(const std::vector<std::string>& list1, const std::vector<std::string>& list2, std::vector<std::string>& overlapping_ids) {
    std::map<std::string, std::pair<int32_t, int32_t> > id2cnt;
    for(size_t i = 0; i < list1.size(); ++i) {
        if ( ++id2cnt[list1[i]].first > 1 ) {
            error("Duplicate ID %s found in list1", list1[i].c_str());
        }
    }
    for(size_t i = 0; i < list2.size(); ++i) {
        if ( ++id2cnt[list2[i]].second > 1 ) {
            error("Duplicate ID %s found in list2", list2[i].c_str());
        }
    }
    overlapping_ids.clear();
    for(std::map<std::string, std::pair<int32_t, int32_t> >::iterator it = id2cnt.begin(); it != id2cnt.end(); ++it) {
        if ( it->second.first > 0 && it->second.second > 0 ) {
            overlapping_ids.push_back(it->first);
        }
    }
    return (int32_t)overlapping_ids.size();
}

// index overlapping IDs from list1 and list2, and store sorted indices in index1 and index2
int32_t index_overlapping_ids(const std::vector<std::string>& list1, const std::vector<std::string>& list2, std::vector<int32_t>& index1, std::vector<int32_t>& index2) {
    std::map<std::string, int32_t> id2idx1;
    std::map<std::string, int32_t> id2idx2;
    for(size_t i = 0; i < list1.size(); ++i) {
        if ( id2idx1.find(list1[i]) == id2idx1.end() ) {
            id2idx1[list1[i]] = i;
        }
        else {
            error("Duplicate ID %s found in list1", list1[i].c_str());
        }
    }
    for(size_t i = 0; i < list2.size(); ++i) {
        if ( id2idx2.find(list2[i]) == id2idx2.end() ) {
            id2idx2[list2[i]] = i;
        }
        else {
            error("Duplicate ID %s found in list2", list2[i].c_str());
        }
    }
    index1.clear();
    index2.clear();
    for(std::map<std::string, int32_t>::iterator it = id2idx1.begin(); it != id2idx1.end(); ++it) {
        if ( id2idx2.find(it->first) != id2idx2.end() ) {
            index1.push_back(it->second);
            index2.push_back(id2idx2[it->first]);
        }
    }
    return (int32_t) index1.size();
}

void fill_list_indices(const std::vector<std::string>& list, const std::map<std::string, int32_t>& list2idx, std::vector<int32_t>& indices)
{
    indices.clear();
    indices.resize(list.size(), -1);
    for (size_t i = 0; i < list.size(); ++i)
    {
        auto it = list2idx.find(list[i]);
        if (it == list2idx.end()) {
            error("Sample ID %s not found in the list provided", list[i].c_str());
        }
        indices[i] = it->second;
    }
}

int32_t rebuild_id2index_map(
    const std::vector<std::string>& ids,
    std::map<std::string, int32_t>& id2idx_map
)
{
    id2idx_map.clear();
    for (size_t i = 0; i < ids.size(); ++i)
    {
        const std::string& id = ids[i];
        if (id2idx_map.find(id) != id2idx_map.end()) {
            error("Duplicate ID %s found when rebuilding ID to index map", id.c_str());
        }
        id2idx_map[id] = i;
    }
    return (int32_t)id2idx_map.size();
}
