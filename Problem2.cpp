#include "basicDS.h"
#include <queue>
#include <algorithm>
#include <set> 
#include <map> 
#define pii pair<int, int> 




/* You can add more functions or variables in each class. 
   But you "Shall Not" delete any functions or variables that TAs defined. */

struct edge{
	int u, v, w; 
	edge() {}
	edge(int u, int v, int w) : u(u), v(v), w(w) {}
	bool operator<(const edge& b) const {
		return w < b.w; 
	}
}; 

struct distance_graph {
	vector<edge> E; 
}; 

struct adj_graph {
	vector<vector<pii>> adj; 
	adj_graph() {}
	adj_graph(Graph& G, int t) {
		adj.resize(G.V.size()+1); 
		for(auto& e : G.E) {
			if(e.b < t) continue; 
			adj[e.vertex[0]].emplace_back(e.vertex[1], e.ce); 
			adj[e.vertex[1]].emplace_back(e.vertex[0], e.ce); 
		}
	}
}; 

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

struct request {
	int s, t, denied_t, id; 
	Set D; 
	bool status; 
	Tree multicastTree; 
	request() {}
	request(int s, int t, int id, Set D, Tree T, bool status) : s(s), t(t), id(id), D(D), multicastTree(T), status(status), denied_t(0) {};
}; 

class Problem2 {
	// vector<vector<pii>> adj; 
public:

	Problem2(Graph G);  //constructor
	~Problem2();        //destructor
	bool insert(int id, int s, Set D, int t, Graph &G, Tree &MTid);
	void stop(int id, Graph &G, Forest &MTidForest);
	void rearrange(Graph &G, Forest &MTidForest);
	void dijkstra(int s, vector<int>& d, vector<int>& p, vector<vector<pii>>& adj); 
	void dfs(int parent, int u, vector<vector<pii>>& adj, vector<vector<int>>& p_v, set<int>& V, Graph& G, Tree& T, int t); 
	void update_time(); 
	map<int, request> rqs; 
};

Problem2::Problem2(Graph G) {
	/* Write your code here. */
	// adj.resize(G.V.size()+1); 
	// for(auto& e : G.E) {
	// 	adj[e.vertex[0]].emplace_back(e.vertex[1], e.ce); 
	// 	adj[e.vertex[1]].emplace_back(e.vertex[0], e.ce); 
	// }
}

Problem2::~Problem2() {
	/* Write your code here. */


}

void Problem2::update_time() {
	for(auto& rq : rqs) {
		if(rq.second.status == false) rq.second.denied_t++; 
	}
}

void Problem2::dijkstra(int s, vector<int>& d, vector<int>& p, vector<vector<pii>>& adj) {
	priority_queue<pii, vector<pii>, greater<pii>> pq;
	d[s] = 0; 
	p[s] = 0; 
	pq.emplace(d[s], s); 
	while(!pq.empty()) {
		
		auto tmp = pq.top(); 
		int disu = tmp.first; 
		int u = tmp.second; 
		pq.pop(); 
		if(d[u] < disu) continue;
		for(auto& nv : adj[u]) {
			
			int v = nv.first; 
			int w = nv.second;
			
			if(d[v] > w + disu) {
				
				d[v] = w + disu; 
				p[v] = u; 
				pq.emplace(d[v], v); 
			}
		}
	}

}

bool check(const graphEdge& a, const graphEdge& b) {
	if(a.vertex[0] != b.vertex[0]) return a.vertex[0] < b.vertex[0]; 
	return a.vertex[1] < b.vertex[1]; 
}

bool equal(const graphEdge& a, const graphEdge& b) {
	return (a.vertex[0] == b.vertex[0] && a.vertex[1] == b.vertex[1]); 
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

int find(graphEdge tar, vector<graphEdge>& a) {
	int pos; 
	pos = search(0, a.size()-1, tar, a).second;
	if(equal(tar, a[pos])) return pos;
	swap(tar.vertex[0], tar.vertex[1]); 
	pos = search(0, a.size()-1, tar, a).second; 
	return pos; 
}

void Problem2::dfs(int parent, int u, vector<vector<pii>>& adj, vector<vector<int>>& p_v, set<int>& V, Graph& G, Tree& T, int t) {
	
	for(auto& nv : adj[u]) {
		int v = nv.first;
		
		if(v != parent) {
			vector<int> path; 
			int cur = v; 
			int tmp_v = v; 
			while(v != 0) {
				path.push_back(v); 
				v = p_v[u][v]; 
			}
			int cnt = 0;
			int first = 0, last = 0;
			for(auto i : path) {
				if(V.find(i) != V.end()) {
					if(!cnt) first = i; 
					cnt++; 
					last = i; 
				}
 			}

			
			if(cnt < 2) {
				for(int i=0;i<path.size()-1;i++) {
					graphEdge ne;
					ne.vertex[0] = path[i]; 
					ne.vertex[1] = path[i+1]; 
					int pos = find(ne, G.E);
					V.insert(path[i]); 
					V.insert(path[i+1]); 
					treeEdge e; 
					e.vertex[0] = G.E[pos].vertex[0]; 
					e.vertex[1] = G.E[pos].vertex[1]; 
					T.E.push_back(e); 
					T.ct += G.E[pos].ce; 
					G.E[pos].b -= t; 
				}
			} else {
				bool flag = false; 
				for(int i=0;i<path.size()-1;i++) {
					if(path[i] == last) flag = false; 
					if(path[i] == first || flag) {
						flag = true; 
						continue;
					} 
					graphEdge ne;
					ne.vertex[0] = path[i]; 
					ne.vertex[1] = path[i+1]; 
					int pos = find(ne, G.E);
					V.insert(path[i]); 
					V.insert(path[i+1]); 
					treeEdge e; 
					e.vertex[0] = G.E[pos].vertex[0]; 
					e.vertex[1] = G.E[pos].vertex[1]; 
					T.E.push_back(e); 
					T.ct += G.E[pos].ce; 
					G.E[pos].b -= t; 
				}
			}
			dfs(u, tmp_v, adj, p_v, V, G, T, t); 
		}
	}
}

bool Problem2::insert(int id, int s, Set D, int t, Graph &G, Tree &MTid) {
	/* Store your output graph and multicast tree into G and MTid */
	if(rqs.find(id) == rqs.end())update_time(); 
	/* Write your code here. */
	// construct GL
	MTid.E.clear(); 
	MTid.V.clear();
	MTid.id = id; 
	MTid.s = s;
	MTid.ct = 0; 
	adj_graph GP(G, t);
	vector<vector<int>> p_v(G.V.size()+1); 
	distance_graph GL; 
	for(auto u : D.destinationVertices) {
		vector<int> p(G.V.size()+1, 0); 
		vector<int> d(G.V.size()+1, INT32_MAX); 
		dijkstra(u, d, p, GP.adj); 
		p_v[u] = p; 
		for(auto& v : D.destinationVertices) {
			if(v != u) {
				if(d[v] == INT32_MAX) {
					rqs.emplace(id, request(s, t, id, D, MTid, false)); 
					return false; 
				}
				GL.E.push_back(edge(u, v, d[v])); 
			}
		}	
	}
	
	// construct TL 
	sort(GL.E.begin(), GL.E.end()); 
	DisjointSet a(G.V.size()); 
	adj_graph TL; 
	TL.adj.resize(G.V.size()+1); 
	for(auto& e : GL.E) {
		if(a.Same(e.u, e.v)) continue;
		a.Union(e.u, e.v); 
		TL.adj[e.u].push_back({e.v, e.w}); 
		TL.adj[e.v].push_back({e.u, e.w}); 	
		
	}
	
	sort(G.E.begin(), G.E.end(), check); 
	set<int> V; 
	dfs(0, D.destinationVertices[0], TL.adj, p_v, V, G, MTid, t); 
	for(auto& u : V) MTid.V.push_back(u); 
	rqs.emplace(id, request(s, t, id, D, MTid, true)); 
	/* You should return true or false according the insertion result */
	return true;
}

void Problem2::stop(int id, Graph &G, Forest &MTidForest) {
	/* Store your output graph and multicast tree forest into G and MTidForest
	   Note: Please "only" include mutlicast trees that you added nodes in MTidForest. */
	update_time(); 
	/* Write your code here. */
	MTidForest.size = 0; 
	MTidForest.trees.clear(); 
	//release bandwidth
	auto& tar = rqs[id]; 
	for(auto& e : tar.multicastTree.E) {
		graphEdge tmp;
		tmp.vertex[0] = e.vertex[0]; 
		tmp.vertex[1] = e.vertex[1]; 
		int pos = find(tmp, G.E);
		G.E[pos].b += tar.t;  
	}
	rqs.erase(id); 

	for(auto& rq : rqs) {
		if(rq.second.status) continue;
		Tree& mt = rq.second.multicastTree; 
		if(insert(mt.id, rq.second.s, rq.second.D, rq.second.t, G, mt)) {
			rq.second.status = true; 
			MTidForest.trees.push_back(mt); 
			MTidForest.size++; 
		}
	}
	
	return;
}

bool com(const Tree& a, const Tree& b) {return a.id < b.id;}

void Problem2::rearrange(Graph &G, Forest &MTidForest) {
	/* Store your output graph and multicast tree forest into G and MTidForest
	   Note: Please include "all" active mutlicast trees in MTidForest. */

	/* Write your code here. */
	MTidForest.size = 0;
	MTidForest.trees.clear(); 
	// first removed all running multicast tree
	for(auto& e : G.E) {
		e.b = e.be; 
	}
	for(auto& rq : rqs) {
		if(rq.second.status == false) continue; 
		Tree& mt = rq.second.multicastTree; 
		if(insert(rq.first, rq.second.s, rq.second.D, rq.second.t, G, mt)) {
			rq.second.status = true; 
			MTidForest.trees.push_back(mt); 
			MTidForest.size++; 
		} else {
			rq.second.status = false; 
		}
	}
	//try to insert new multicast tree
	auto compare = [](request a, request b) {return a.denied_t > b.denied_t;}; 

	set<request, decltype(compare)> unsatisfied(compare); 
	for(auto& rq : rqs) {
		if(rq.second.status == false) unsatisfied.insert(rq.second); 
		
	}

	for(auto& rq : unsatisfied) {
		Tree tmp;
		if(insert(rq.id, rq.s, rq.D, rq.t, G, tmp)) {
			rqs[rq.id].denied_t = 0; 
			rqs[rq.id].multicastTree = tmp; 
			rqs[rq.id].status = true; 
			MTidForest.trees.push_back(tmp); 
			MTidForest.size++; 
		}
	}
	sort(MTidForest.trees.begin(), MTidForest.trees.end(), com); 	
	update_time(); 
	return;
}
