#ifndef TURBO_Utility_DSUClusteredSet
#define TURBO_Utility_DSUClusteredSet

#include <vector>
#include <functional>
#include <unordered_set>
#include <map>

template<class T>
class DSUClusteredSet {
private:
	//Size of set
	size_t _size;

	//Storage for set elements
	std::vector<T> _elements;

	//DSU implementation inside
	std::vector<int> _rank;
	std::vector<int> _parent;
public:
	void InitSet(std::vector<T>& elements) {
		_elements = elements;
		_size = elements.size();
		_rank.resize(_size);
		_parent.resize(_size);
		for (int i = 0; i < _size; i++) {
			MakeSet(i);
		};
	};

	void MakeSet(int x) {
		_rank[x] = 0;
		_parent[x] = x;
	};

	int Find(int x) {
		if (x != _parent[x]) {
			_parent[x] = Find(_parent[x]);
		};
		return _parent[x];
	};

	void Union(int x, int y) {
		int xRoot = Find(x);
		int yRoot = Find(y);
		if (xRoot == yRoot) return;

		// x and y are not already in same set. Merge them.
		if (_rank[xRoot] < _rank[yRoot]) {
			_parent[xRoot] = yRoot;
		} else if (_rank[xRoot] > _rank[yRoot]) {
			_parent[yRoot] = xRoot;
		} else {
			_parent[yRoot] = xRoot;
			_rank[xRoot] = _rank[xRoot] + 1;
		};
	};

	//Where given predicate equals true
	std::vector<std::vector<T*> > GetClusters( std::function<bool(T)> P ) {
		//Determine number of different root elements
		std::unordered_set<int> rootsSet;
		for (int i = 0; i<_size; i++) {
			if (!P(_elements[i])) continue;
			rootsSet.insert(Find(i));
		};

		//For each cluster assemble elements
		size_t nClusters = rootsSet.size();	
		std::vector<int> roots;	
		std::map<int, int> rootToClusterIndex;
		std::vector< std::vector<T*> > clusterIndexToElements;
		int clusterIndex = 0;
		for (int root : rootsSet) {
			roots.push_back(root);
			std::vector<T*> clusterElements;
			clusterIndexToElements.push_back(  clusterElements );
			rootToClusterIndex[root] = clusterIndex++;
		};
		
		//Collect elements
		for (int i = 0; i<_size; i++) {
			int root = Find(i);
			if (rootsSet.find(root) != std::end(rootsSet)) {
				int cI = rootToClusterIndex[root];
				clusterIndexToElements[cI].push_back(&_elements[i]);
			};
		};

		//Output clusters
		return clusterIndexToElements;
	};
};

#endif
