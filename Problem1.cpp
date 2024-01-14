#include "basicDS.h"
#include <algorithm>
#include <set>
#include <map> 

/* You can add more functions or variables in each class. 
   But you "Shall Not" delete any functions or variables that TAs defined. */
using pii = pair<int, int>; 
class DisjointSet {
    private : 
    vector<int> parent;
    int size; 
    public : 
    int find_root(int x) {
        if(x == parent[x]) return x;
        int root = find_root(parent[x]);
        return parent[x] = root;
    }
    DisjointSet(int n) : parent(n+1), size(n) {init();}
    void init() {
        for(int i=0;i<parent.size();i++) {
            parent[i] = i;
        }
    }
    bool Same(int a, int b) {return find_root(a) == find_root(b);}
    void Union(int a, int b) {
        int roota = find_root(a);
        int rootb = find_root(b);
        if(roota != rootb) size--;
        parent[roota] = find_root(rootb);
    }
    int get_size() {return size;}
};

bool compare(const graphEdge& a,const graphEdge& b) {
	return a.ce < b.ce; 
}

bool T_compare(const Tree& a, const Tree& b) {
	return a.id < b.id; 
}

struct request {
	int s, t; 
	Set D; 
	bool status; 
	Tree multicastTree; 
	request() {}
	request(int s, int t, Set D, Tree T) : s(s), t(t), D(D), multicastTree(T) {status = true;};
}; 

class Problem1 {
public:

	Problem1(Graph G);  //constructor
	~Problem1();        //destructor
	void insert(int id, int s, Set D, int t, Graph &G, Tree &MTid);
	void stop(int id, Graph &G, Forest &MTidForest);
	void rearrange(Graph &G, Forest &MTidForest);
	// vector<graphEdge> G_W, G_ID; 
	map<int, request> rqs; 
};

Problem1::Problem1(Graph G) {
	/* Write your code here. */
	// G_W = G.E; 
	// G_ID = G.E; 
	// sort(G_W.begin(), G_W.end(), compare); 
	// sort(G_ID.begin(), G_ID.end(), check); 
}

Problem1::~Problem1() {
	/* Write your code here. */

}

void Problem1::insert(int id, int s, Set D, int t, Graph &G, Tree &MTid) {
	/* Store your output graph and multicast tree into G and MTid */
	
	/* Write your code here. */
	sort(G.E.begin(), G.E.end(), compare);
	DisjointSet a(G.V.size());
	vector<bool> take(G.E.size(), false);
	for(int i=0;i<G.E.size();i++) {
		if(G.E[i].b < t || a.Same(G.E[i].vertex[0], G.E[i].vertex[1])) continue;
		a.Union(G.E[i].vertex[0], G.E[i].vertex[1]); 
		take[i] = true;
	}
	MTid.E.clear();
	MTid.V.clear();
	set<int> tmp_v;
	tmp_v.insert(s); 
	MTid.ct = 0;
	for(int i=0;i<G.E.size();i++) {
		if(!take[i]) continue;
		if(a.Same(s, G.E[i].vertex[0])) {
			G.E[i].b -= t;
			treeEdge ne;
			ne.vertex[0] = G.E[i].vertex[0];
			ne.vertex[1] = G.E[i].vertex[1];
			MTid.E.push_back(ne);
			tmp_v.insert(ne.vertex[0]);
			tmp_v.insert(ne.vertex[1]);
			MTid.ct += G.E[i].ce; 
		}
	}
	for(auto& u : tmp_v) MTid.V.push_back(u); 
	MTid.id = id;
	MTid.s = s; 
	// runingTrees.trees.push_back(MTid);
	// runingTrees.size++;
	rqs.emplace(id, request(s, t, D, MTid)); 
	return;
}

bool check(const graphEdge& a, const graphEdge& b) {
	if(a.vertex[0] != b.vertex[0]) return a.vertex[0] < b.vertex[0]; 
	return a.vertex[1] < b.vertex[1]; 
}

pii search(int L, int R, graphEdge tar, vector<graphEdge>& a) {
	if(check(a[L], tar) == false) return {L-1, L}; 
	while(R-L>1) {
		int m = L + (R - L)/2; 
		if(check(a[m], tar)) 
			L = m; 
		else 
			R = m;
 	}
	return {L, R}; 
}


void Problem1::stop(int id, Graph &G, Forest &MTidForest) {
	/* Store your output graph and multicast tree forest into G and MTidForest
	   Note: Please "only" include mutlicast trees that you added nodes in MTidForest. */
	
	/* Write your code here. */
	//clear forest
	MTidForest.size = 0;
	MTidForest.trees.clear(); 
	//sort graph edge for binary search
	sort(G.E.begin(), G.E.end(), check);
	
	// release the bandwidth
	auto& tar = rqs[id]; 
	for(auto& e : tar.multicastTree.E) {
		graphEdge tmp;
		tmp.vertex[0] = e.vertex[0]; 
		tmp.vertex[1] = e.vertex[1]; 
		int pos = search(0, G.E.size()-1, tmp, G.E).second;
		G.E[pos].b = min(G.E[pos].be, tar.t + G.E[pos].b);  
	}
	rqs.erase(id); 
	
	for(auto& rt : rqs) {
		vector<graphEdge> E_inw; 
		for(auto& e : G.E) E_inw.push_back(e); 
		sort(E_inw.begin(), E_inw.end(), compare); 
		auto& mt = rt.second.multicastTree; 
		if(mt.V.size() == G.V.size()) continue; //already connect all vertex.
		DisjointSet a(G.V.size()); 
		vector<bool> take(G.E.size(), false);
		vector<bool> original(G.E.size(), false); 

		for(auto& e : mt.E) {
			a.Union(e.vertex[0], e.vertex[1]); 
			graphEdge tmp;
			tmp.vertex[0] = e.vertex[0];
			tmp.vertex[1] = e.vertex[1];
			int pos = search(0, G.E.size()-1, tmp, G.E).second; 
			take[pos] = true; 
			original[pos] = true; 
		}
		int traffic = rt.second.t; 
		bool update = false; 
		// running kruskal 
		for(int i=0;i<E_inw.size();i++) {
			if(E_inw[i].b < traffic || a.Same(E_inw[i].vertex[0], E_inw[i].vertex[1])) continue;
			a.Union(E_inw[i].vertex[0], E_inw[i].vertex[1]); 
			int pos = search(0, E_inw.size()-1, E_inw[i], G.E).second; 
			take[pos] = true; 
		} 
		// add new edge and vertex to tree 
		set<int> tmp_v; 
		for(auto& u : mt.V) tmp_v.insert(u); 
		for(int i=0;i<G.E.size();i++) {
			if(!take[i] || !a.Same(mt.s, G.E[i].vertex[0]) || original[i]) continue;
			update = true; 
			G.E[i].b -= traffic; 
			treeEdge ne;
			ne.vertex[0] = G.E[i].vertex[0];
			ne.vertex[1] = G.E[i].vertex[1]; 
			mt.E.push_back(ne); 
			mt.ct += G.E[i].ce; 
			tmp_v.insert(ne.vertex[0]);
			tmp_v.insert(ne.vertex[1]); 
		}
		mt.V.clear(); 
		for(auto& u : tmp_v) mt.V.push_back(u); 
		if(update) {
			MTidForest.trees.push_back(mt); 
			MTidForest.size++;
		}
	}

	return;
}

void Problem1::rearrange(Graph &G, Forest &MTidForest) {
	/* Store your output graph and multicast tree forest into G and MTidForest
	   Note: Please include "all" active mutlicast trees in MTidForest. */

	/* Write your code here. */
	MTidForest.size = 0; 
	MTidForest.trees.clear(); 
	for(auto& e : G.E) {
		e.b = e.be; 
	}
	for(auto& rt : rqs) {
		
		Tree& mt = rt.second.multicastTree; 
		insert(rt.first, rt.second.s, rt.second.D, rt.second.t, G, mt);
		MTidForest.trees.push_back(mt); 
		MTidForest.size++; 
	}
	return;
}
