using namespace std;

int *loc, *tmploc, cvCnt = 0;
class edge
{
public:
	int s;
	int t;
};

class idW
{
public:
	int id;
	int weight;
};

inline int max3(int a, int b, int c)
{
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

bool cmp(const pair<int, int> &a, const pair<int, int> &b)
{
	if (a.second == b.second)
		return a.first < b.first;
	else
		return a.second > b.second;
}
bool IGCmp(const idW &a, const idW &b)
{
	return a.weight == b.weight ? (a.id > b.id) : (a.weight > b.weight);
}

class Clique
{
public:
	int *lab, *cdv, *adjv, *ns, **degS, **subS, *ver, *dTMP, maxK, curK;
	long long cliqueNumber;
	Clique(int n, int e, int k);
	~Clique();
};

Clique::Clique(int n, int e, int k)
{
	lab = new int[n]();
	dTMP = new int[n];
	cdv = new int[n + 1];
	adjv = new int[e * 2];
	ns = new int[k + 1];
	ver = new int[k + 1];
	degS = new int *[k + 1], subS = new int *[k + 1];
	for (int i = 1; i <= k; i++)
	{
		degS[i] = new int[n];
		subS[i] = new int[n];
	}
	maxK = k;
	cliqueNumber = 0;
}
Clique::~Clique(void)
{
	// printf("delete clique\n");
	if (lab != NULL)
		delete[] lab;
	if (dTMP != NULL)
		delete[] dTMP;
	if (cdv != NULL)
		delete[] cdv;
	if (adjv != NULL)
		delete[] adjv;
	if (ns != NULL)
		delete[] ns;
	if (ver != NULL)
		delete[] ver;

	if (degS != NULL)
		for (int i = 1; i <= maxK; i++)
			if (degS[i] != NULL)
				delete[] degS[i];
	delete[] degS;

	if (subS != NULL)
		for (int i = 1; i <= maxK; i++)
			if (subS[i] != NULL)
				delete[] subS[i];
	delete[] subS;
}

class Graph
{
public:
	Graph();
	Graph(const Graph &obj);
	~Graph();
	void readedgelist(char *edgelist);
	void coreDecomposition();
	void mkGraph();
	int outLargeClique();
	bool isEdge(int, int);
	Graph &mksub(int, ...);
	// Graph* mksubMark(int*, int, int*);
	int color(int *, int);
	int Greedy_SD(int *color);
	void kClique(int, long long *, long long *);
	void kCliqueCount(int, long long *, int *, int *, int *, int *, int *, int **, int **, long long *, int);

	void kCliqueNew(int k, long long *tol, long long *cnt, int *subg, int size);
	void kCliqueCountNew(int l, long long *tol, long long *cnt);

	int deleteNodes(int *delArray, int size);
	void deleteEdge(int s, int t);
	void insertEdge(int s, int t);
	// bool isEdge(int, int);

	int n;
	int e;
	int maxDeg;
	edge *edges;

	int *deg;
	int *cd;
	int *adj;
	int *coreRank; // increasing core number order
	int *coreNum;  // coreNum[i] is the core number of node i.
	int *bin;

	Clique *clique;
};

Graph::Graph(void)
{
	edges = NULL;
	deg = cd = adj = coreRank = coreNum = bin = NULL;
	clique = NULL;
}
Graph::~Graph(void)
{
	if (edges != NULL)
		delete[] edges;
	if (deg != NULL)
		delete[] deg;
	if (cd != NULL)
		delete[] cd;
	if (adj != NULL)
		delete[] adj;
	if (coreRank != NULL)
		delete[] coreRank;
	if (coreNum != NULL)
		delete[] coreNum;
	if (bin != NULL)
		delete[] bin;
	if (clique != NULL)
		delete clique;
}
Graph::Graph(const Graph &obj)
{
	n = obj.n, e = obj.e, maxDeg = obj.maxDeg, edges = obj.edges;
	edges = obj.edges;

	if (deg != NULL)
		delete[] deg;
	if (obj.deg != NULL)
		deg = new int[n], memcpy(deg, obj.deg, n * sizeof(int));

	if (cd != NULL)
		delete[] cd;
	if (obj.cd != NULL)
		cd = new int[n + 1], memcpy(cd, obj.cd, (n + 1) * sizeof(int));

	if (adj != NULL)
		delete[] adj;
	if (obj.adj != NULL)
		adj = new int[2 * e], memcpy(adj, obj.adj, 2 * e * sizeof(int));

	if (coreRank != NULL)
		delete[] coreRank;
	if (obj.coreRank != NULL)
		coreRank = new int[n], memcpy(coreRank, obj.coreRank, n * sizeof(int));

	if (coreNum != NULL)
		delete[] coreNum;
	if (obj.coreNum != NULL)
		coreNum = new int[n], memcpy(coreNum, obj.coreNum, n * sizeof(int));

	if (bin != NULL)
		delete[] bin;
	if (obj.bin != NULL)
		bin = new int[maxDeg + 2], memcpy(bin, obj.bin, (maxDeg + 2) * sizeof(int));
}

class FastRead
{
	char *inArr;
	char *fileName;
	int index;

public:
	FastRead(char *fileName) : fileName(fileName), index(0){};

	void FileToArr();
	inline int read();
	~FastRead();
};

void FastRead::FileToArr()
{
	FILE *pf = fopen(fileName, "r");
	fseek(pf, 0, SEEK_END);
	long long lSize = ftell(pf);

	inArr = new char[lSize + 1];
	// printf("lSize = %lld\n", lSize);

	rewind(pf);
	fread(inArr, sizeof(char), lSize, pf);
	fclose(pf);
	inArr[lSize] = '\0';
}

inline int FastRead::read()
{
	int s = 0, w = 1;
	char ch = inArr[index++];
	while (ch < '0' || ch > '9')
	{
		if (ch == '-')
			w = -1;
		ch = inArr[index++];
	}
	while (ch >= '0' && ch <= '9')
		s = s * 10 + ch - '0', ch = inArr[index++];
	return s * w;
}

FastRead::~FastRead()
{
	delete[] inArr;
}

void Graph::readedgelist(char *edgelist)
{
	FastRead fr(edgelist);
	fr.FileToArr();

	printf("input finished.\n");

	n = fr.read();
	e = fr.read();

	int edgeN = e;

	edges = new edge[e];
	e = 0;

	while (e < edgeN)
	{
		edges[e].s = fr.read();
		edges[e].t = fr.read();
		e++;
	}
}

void Graph::mkGraph()
{
	deg = new int[n]();
	cd = new int[n + 1];
	adj = new int[2 * e];
	maxDeg = 0;
	for (int i = 0; i < e; i++)
	{
		deg[edges[i].s]++;
		deg[edges[i].t]++;
		maxDeg = max3(maxDeg, deg[edges[i].s], deg[edges[i].t]);
	}
	cd[0] = 0;
	for (int i = 1; i < n + 1; i++)
	{
		cd[i] = cd[i - 1] + deg[i - 1];
		deg[i - 1] = 0;
	}

	for (int i = 0; i < e; i++)
	{
		adj[cd[edges[i].s] + deg[edges[i].s]++] = edges[i].t;
		adj[cd[edges[i].t] + deg[edges[i].t]++] = edges[i].s;
	}
	// delete [] edges;

	for (int i = 0; i < n; i++)
		sort(adj + cd[i], adj + cd[i] + deg[i]);
}

bool Graph::isEdge(int a, int b)
{
	if (deg[a] > deg[b])
		a = a ^ b, b = a ^ b, a = a ^ b;
	for (int i = cd[a]; i < cd[a] + deg[a]; i++)
		if (adj[i] == b)
			return true;
	return false;
}
int Graph::outLargeClique()
{
	int CSize = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		int id = coreRank[i];
		if (coreNum[id] >= CSize)
		{
			pair<int, int> *SCore = new pair<int, int>[deg[id]];
			int cnt = 0, ind = 0;
			for (int j = cd[id]; j < cd[id] + deg[id]; j++)
				if (coreNum[adj[j]] >= CSize)
				{
					SCore[cnt].first = adj[j];
					SCore[cnt].second = coreNum[adj[j]];
					cnt++;
				}
			sort(SCore, SCore + cnt, cmp);
			int *C = new int[deg[id]];
			for (int j = 0; j < cnt; j++)
			{
				int flag = 1;
				for (int k = 0; k < ind; k++)
				{
					if (isEdge(SCore[j].first, C[k]) == false)
					{
						flag = 0;
						break;
					}
				}
				if (flag)
					C[ind++] = SCore[j].first;
			}
			ind++; // node "id" ?
			if (ind > CSize)
				CSize = ind;
			delete[] SCore;
		}
	}
	return CSize;
}

void Graph::coreDecomposition()
{
	bin = new int[maxDeg + 2]();

	for (int i = 0; i < n; i++)
		bin[deg[i]]++;

	int lastBin = bin[0], nowBin;
	bin[0] = 0;
	for (int i = 1; i <= maxDeg; i++)
	{
		nowBin = lastBin + bin[i - 1];
		lastBin = bin[i];
		bin[i] = nowBin;
	}
	int *vert = new int[n](), *pos = new int[n](), *tmpDeg = new int[n]();
	for (int i = 0; i < n; i++)
	{
		pos[i] = bin[deg[i]];

		vert[bin[deg[i]]++] = i;
		tmpDeg[i] = deg[i];
	}

	bin[0] = 0;
	for (int i = maxDeg; i >= 1; i--)
	{
		bin[i] = bin[i - 1];
	}

	int *cNum = new int[n];
	for (int i = 0; i < n; i++)
	{
		int id = vert[i], nbr, binFrontId;
		cNum[id] = tmpDeg[id];
		for (int i = cd[id]; i < cd[id] + deg[id]; i++)
		{
			nbr = adj[i];

			if (tmpDeg[nbr] > tmpDeg[id])
			{
				binFrontId = vert[bin[tmpDeg[nbr]]];
				if (binFrontId != nbr)
				{

					pos[binFrontId] = pos[nbr];
					pos[nbr] = bin[tmpDeg[nbr]];
					vert[bin[tmpDeg[nbr]]] = nbr;
					vert[pos[binFrontId]] = binFrontId;
				}
				bin[tmpDeg[nbr]]++;
				tmpDeg[nbr]--;
			}
		}
	}
	coreNum = cNum;
	coreRank = vert;
	delete[] tmpDeg;
	delete[] pos;
}

Graph &Graph::mksub(int argCnt, ...) //(int *nodes, int NodeNum)//
{
	Graph *sg = new Graph;
	int *Mark = NULL;
	va_list ap;
	va_start(ap, argCnt);
	// int *sg = va_arg(ap,int *);

	int *nodes = va_arg(ap, int *);
	int NodeNum = va_arg(ap, int);

	if (argCnt >= 3)
		Mark = va_arg(ap, int *);

	va_end(ap);

	if (NodeNum < n / 1000)
	{
		sg->n = NodeNum, sg->e = 0;
		int *lab = new int[n], cnt = 0;
		for (int i = 0; i < sg->n; i++)
		{
			lab[nodes[i]] = cnt++;
		}

		if (argCnt >= 3)
		{
			cnt = 0;
			for (int i = 0; i < sg->n; i++)
				Mark[cnt++] = nodes[i];
		}

		for (int i = 0; i < NodeNum; i++)
		{
			for (int j = i + 1; j < NodeNum; j++)
			{
				for (int p = cd[nodes[i]]; p < cd[nodes[i] + 1]; p++)
				{
					if (adj[p] == nodes[j])
					{
						sg->e++;
					}
				}
			}
		}
		sg->edges = new edge[sg->e];
		sg->e = 0;

		for (int i = 0; i < NodeNum; i++)
		{
			for (int j = i + 1; j < NodeNum; j++)
			{
				for (int p = cd[nodes[i]]; p < cd[nodes[i] + 1]; p++)
				{
					if (adj[p] == nodes[j])
					{
						sg->edges[sg->e].s = lab[nodes[i]];
						sg->edges[sg->e].t = lab[nodes[j]];
						sg->e++;
					}
				}
			}
		}
		delete[] lab;
		sg->mkGraph();
		return *sg;
	}

	sg->n = NodeNum, sg->e = 0;
	int *newFg = new int[n]();

	for (int i = 0; i < NodeNum; i++)
		newFg[nodes[i]] = 1;

	for (int i = 0; i < e; i++)
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1)
			sg->e++;

	sg->edges = new edge[sg->e];

	sg->e = 0;
	int *lab = new int[n], cnt = 0;
	for (int i = 0; i < sg->n; i++)
	{
		lab[nodes[i]] = cnt++;
	}

	if (argCnt >= 3)
	{
		cnt = 0;
		for (int i = 0; i < sg->n; i++)
			Mark[cnt++] = nodes[i];
	}

	for (int i = 0; i < e; i++)
	{
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1)
		{
			sg->edges[sg->e].s = lab[edges[i].s];
			sg->edges[sg->e].t = lab[edges[i].t];
			sg->e++;
		}
	}

	delete[] newFg;
	delete[] lab;
	sg->mkGraph();
	return *sg;
}

class SDNode
{
public:
	int adjColorsNum;
	int adjUncolored;

	SDNode(){};
	SDNode(int c, int u) : adjColorsNum(c), adjUncolored(u){};

	~SDNode(){};

	bool operator<(SDNode &sdn)
	{
		if (this->adjColorsNum > sdn.adjColorsNum) // here using ">" rather than "<" cause we only build Min-heap in heapLLU.h
			return true;
		if (this->adjColorsNum == sdn.adjColorsNum && this->adjUncolored > sdn.adjUncolored)
			return true;
		return false;
	}
	bool operator>(SDNode &sdn)
	{
		return !(*this < sdn);
	}
	bool operator==(const SDNode &sdn)
	{
		return this->adjColorsNum == sdn.adjColorsNum && this->adjUncolored == sdn.adjUncolored;
	}
	SDNode &operator=(SDNode &sdn)
	{
		if (*this == sdn)
			return *this;

		this->adjColorsNum = sdn.adjColorsNum;
		this->adjUncolored = sdn.adjUncolored;
		return *this;
	}
};

int Graph::Greedy_SD(int *color)
{
	memset(color, -1, sizeof(int) * n);
	SDNode *nodes = new SDNode[n];
	for (int i = 0; i < n; i++)
	{
		nodes[i].adjColorsNum = 0;
		nodes[i].adjUncolored = deg[i];
	}
	bheapLLU<SDNode> *SDHeap = mkheapLLU<SDNode>(n, nodes);
	// bool* nodeFlag = new bool[n]();

	int leftNodes = n;
	keyvalueLLU<SDNode> kv, kvNew;

	bool *C = new bool[maxDeg + 1]();

	int nodeID, nbrID, nnbrID;
	int colorNum = 1;

	while (leftNodes > 0)
	{
		kv = popminLLU<SDNode>(SDHeap);
		nodeID = kv.key;

		for (int i = cd[nodeID]; i < cd[nodeID] + deg[nodeID]; i++)
		{
			nbrID = adj[i];
			if (color[nbrID] != -1)
				C[color[nbrID]] = 1;
		}
		for (int i = 0; i < maxDeg + 1; i++)
			if (C[i] == 0)
			{
				color[nodeID] = i;
				colorNum = i > colorNum ? i : colorNum;
				break;
			}

		for (int i = cd[nodeID]; i < cd[nodeID] + deg[nodeID]; i++)
		{
			nbrID = adj[i];
			if (color[nbrID] != -1)
				C[color[nbrID]] = 0;
		}

		// nodeFlag[nodeID] = 1;

		int nbrColorNum;
		for (int i = cd[nodeID]; i < cd[nodeID] + deg[nodeID]; i++)
		{
			nbrID = adj[i];
			nbrColorNum = 0;
			if (color[nbrID] == -1) // uncolored
			{
				for (int j = cd[nbrID]; j < cd[nbrID] + deg[nbrID]; j++)
				{
					nnbrID = adj[j];

					if (color[nnbrID] != -1 && C[color[nnbrID]] == 0)
					{
						C[color[nnbrID]] = 1;
						nbrColorNum++;
					}
				}
				// kvNew.key = nbrID;
				// kvNew.value.adjColorsNum = nbrColorNum;
				// kvNew.value.adjUncolored = kv.value.adjUncolored - 1;

				updateLLU(SDHeap, nbrID, SDNode(nbrColorNum, kv.value.adjUncolored - 1));

				for (int j = cd[nbrID]; j < cd[nbrID] + deg[nbrID]; j++)
				{
					nnbrID = adj[j];
					if (color[nnbrID] != -1)
						C[color[nnbrID]] = 0;
				}
			}
		}

		leftNodes--;
	}
	cout << "------------------------" << endl;

	return colorNum + 1;
}

int Graph::color(int *color, int colorAlgo = 0)
{
	cout << "colorAlgo: " << colorAlgo << endl;
	if (colorAlgo == 3) // Greedy SD coloring algorithm
	{
		return Greedy_SD(color);
	}

	idW *ig = new idW[n];

	if (colorAlgo == 0) // degree order
	{
		for (int i = 0; i < n; i++)
		{
			ig[i].id = i;
			ig[i].weight = deg[i];
		}
	}

	if (colorAlgo == 1) // degeneracy order
	{
		if (coreNum == NULL)
		{
			coreDecomposition();
		}

		for (int i = 0; i < n; i++)
		{
			ig[i].id = i;
			ig[i].weight = coreNum[i];
		}
	}

	if (colorAlgo == 2) // first fit order
	{
		for (int i = 0; i < n; i++)
		{
			ig[i].id = i;
			ig[i].weight = i;
		}
	}

	sort(ig, ig + n, IGCmp);

	memset(color, -1, sizeof(int) * n);
	int *C = new int[(maxDeg + 1)]();

	color[ig[0].id] = 0;
	int colorNum = 1;

	for (int i = 1; i < n; i++)
	{
		int tmpDeg = deg[ig[i].id], tmpid = ig[i].id;
		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 1;
		}
		for (int j = 0; j < maxDeg + 1; j++)
			if (C[j] == 0)
			{
				color[ig[i].id] = j;
				colorNum = j > colorNum ? j : colorNum;
				break;
			}

		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 0;
		}
	}
	delete[] ig;
	delete[] C;
	return colorNum + 1;
}

void Graph::kCliqueCount(int l, long long *tol,
						 int *ver, int *lab, int *cdv, int *adjv, int *ns, int **degS, int **subS, long long *cnt, int K)
{
	int u, v, w, end;
	if (l == 2)
	{
		for (int i = 0; i < ns[2]; i++)
		{ // list all edges
			u = subS[2][i];
			ver[2] = u;
			end = cdv[u] + degS[2][u];
			for (int p = 2; p <= K; p++)
				cnt[ver[p]] += degS[2][u];

			(*tol) += degS[2][u];

			for (int j = cdv[u]; j < end; j++)
			{
				// listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
				cnt[adjv[j]]++;
			}
		}
		return;
	}

	for (int i = 0; i < ns[l]; i++)
	{
		u = subS[l][i];
		ver[l] = u;
		ns[l - 1] = 0;
		end = cdv[u] + degS[l][u];
		for (int j = cdv[u]; j < end; j++) // relabeling nodes and forming U'.
		{
			v = adjv[j];
			if (lab[v] == l)
			{
				lab[v] = l - 1;
				subS[l - 1][ns[l - 1]++] = v;
				degS[l - 1][v] = 0; // new degrees
			}
		}
		for (int j = 0; j < ns[l - 1]; j++) // reodering adjacency list and computing new degrees
		{
			v = subS[l - 1][j];
			end = cdv[v] + degS[l][v];
			for (int k = cdv[v]; k < end; k++)
			{
				w = adjv[k];
				if (lab[w] == l - 1)
					degS[l - 1][v]++;
				else
				{
					adjv[k--] = adjv[--end];
					adjv[end] = w;
				}
			}
		}

		kCliqueCount(l - 1, tol, ver, lab, cdv, adjv, ns, degS, subS, cnt, K);

		for (int j = 0; j < ns[l - 1]; j++)
		{ // restoring labels
			v = subS[l - 1][j];
			lab[v] = l;
		}
	}
}

void Graph::kClique(int k, long long *tol, long long *cnt)
{
	int *d, *sub, *lab, *cdv, *adjv, *ns, **degS, **subS;
	int nsg, maxDv;

	for (int i = 0; i < e; i++)
	{
		if (deg[edges[i].s] > deg[edges[i].t])
		{
			swap(edges[i].s, edges[i].t);
		}
	}

	// mkspecial

	d = new int[n]();
	for (int i = 0; i < e; i++)
		d[edges[i].s]++;

	cdv = new int[n + 1];
	nsg = 0, maxDv = 0, cdv[0] = 0;
	sub = new int[n], lab = new int[n];

	for (int i = 1; i < n + 1; i++)
	{
		cdv[i] = cdv[i - 1] + d[i - 1];
		maxDv = (maxDv > d[i - 1]) ? maxDv : d[i - 1];
		sub[nsg++] = i - 1;
		d[i - 1] = 0;
		lab[i - 1] = k;
	}

	adjv = new int[e];
	for (int i = 0; i < e; i++)
		adjv[cdv[edges[i].s] + d[edges[i].s]++] = edges[i].t;

	ns = new int[k + 1];
	ns[k] = nsg;

	degS = new int *[k + 1], subS = new int *[k + 1];

	for (int i = 2; i < k; i++)
	{
		degS[i] = new int[n];
		subS[i] = new int[maxDv];
	}
	degS[k] = d;
	subS[k] = sub;

	int *ver = new int[k + 1];
	kCliqueCount(k, tol, ver, lab, cdv, adjv, ns, degS, subS, cnt, k);

	delete[] lab, delete[] cdv, delete[] adjv, delete[] ns;

	for (int i = 2; i <= k; i++)
		delete[] degS[i], delete[] subS[i];
	delete[] degS, delete[] subS;
}

void Graph::kCliqueCountNew(int l, long long *tol, long long *cnt)
{
	int u, v, w, end;
	if (l == 1)
	{
		(*tol) += clique->ns[1];
		int k = clique->curK;

		for (int p = 2; p <= k; p++)
			cnt[clique->ver[p]] += clique->ns[1];

		for (int i = 0; i < clique->ns[1]; i++)
		{
			u = clique->subS[1][i];
			cnt[u]++;
			clique->ver[1] = u;
		}
		return;
	}
	if (l == 2)
	{
		int k = clique->curK;
		for (int i = 0; i < clique->ns[2]; i++)
		{ // list all edges
			u = clique->subS[2][i];
			clique->ver[2] = u;
			end = clique->cdv[u] + clique->degS[2][u];
			for (int p = 2; p <= k; p++)
				cnt[clique->ver[p]] += clique->degS[2][u];

			(*tol) += clique->degS[2][u];

			for (int j = clique->cdv[u]; j < end; j++)
			{
				// listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
				cnt[clique->adjv[j]]++;
			}
		}
		return;
	}

	for (int i = 0; i < clique->ns[l]; i++)
	{
		u = clique->subS[l][i];
		clique->ver[l] = u;
		clique->ns[l - 1] = 0;
		end = clique->cdv[u] + clique->degS[l][u];
		for (int j = clique->cdv[u]; j < end; j++) // relabeling nodes and forming U'.
		{

			v = clique->adjv[j];
			if (clique->lab[v] == l)
			{
				clique->lab[v] = l - 1;
				clique->subS[l - 1][clique->ns[l - 1]++] = v;
				clique->degS[l - 1][v] = 0; // new degrees
			}
		}
		for (int j = 0; j < clique->ns[l - 1]; j++) // reodering adjacency list and computing new degrees
		{
			v = clique->subS[l - 1][j];
			end = clique->cdv[v] + clique->degS[l][v];
			for (int k = clique->cdv[v]; k < end; k++)
			{
				w = clique->adjv[k];
				if (clique->lab[w] == l - 1)
					clique->degS[l - 1][v]++;
				else
				{
					clique->adjv[k--] = clique->adjv[--end];
					clique->adjv[end] = w;
				}
			}
		}
		kCliqueCountNew(l - 1, tol, cnt);

		for (int j = 0; j < clique->ns[l - 1]; j++)
		{ // restoring labels
			v = clique->subS[l - 1][j];
			clique->lab[v] = l;
		}
	}
}

// cd, adj, deg, subS(subgraph node set), ns(subgraph node size)
void Graph::kCliqueNew(int k, long long *tol, long long *cnt, int *subg, int size)
{
	int cdTmp = 0, nsg = 0;
	for (int i = 0; i < size; i++)
		clique->lab[subg[i]] = k;
	for (int i = 0; i < size; i++)
	{
		int u = subg[i];
		clique->degS[k][u] = 0;
		for (int j = cd[u]; j < cd[u] + deg[u]; j++)
		{
			int v = adj[j];
			if (clique->lab[v] == k)
				clique->degS[k][u]++;
		}
	}
	for (int i = 0; i < size; i++)
	{
		int u = subg[i];
		clique->dTMP[u] = 0;
		clique->cdv[u] = cd[u];
		clique->subS[k][nsg++] = u;
		for (int j = cd[u]; j < cd[u] + deg[u]; j++)
		{
			int v = adj[j];
			if (clique->lab[v] == k)
			{
				if (clique->degS[k][u] < clique->degS[k][v] || (clique->degS[k][u] == clique->degS[k][v] && u < v))
				{
					clique->adjv[clique->cdv[u] + clique->dTMP[u]++] = v;
				}
			}
		}
	}

	for (int i = 0; i < size; i++)
		clique->degS[k][subg[i]] = clique->dTMP[subg[i]];

	clique->ns[k] = nsg;
	clique->curK = k;
	kCliqueCountNew(k, tol, cnt);
	if (size == n)
		clique->cliqueNumber = *tol;
	for (int i = 0; i < size; i++)
		clique->lab[subg[i]] = 0;
}

int Graph::deleteNodes(int *delArray, int size)
{
	if (size == 1)
	{
		int u = delArray[0];
		for (int j = cd[u]; j < cd[u] + deg[u]; j++)
		{
			int v = adj[j];
			for (int k = cd[v]; k < cd[v] + deg[v]; k++)
			{
				int w = adj[k];
				if (w == u)
				{
					adj[k] = adj[cd[v] + deg[v] - 1];
					adj[cd[v] + deg[v] - 1] = w;
					deg[v]--;
					break;
				}
			}
		}
		deg[u] = 0;
		return 0;
	}

	bool *delMark = new bool[n]();

	for (int i = 0; i < size; i++)
		delMark[delArray[i]] = true;

	int interEdges = 0;
	for (int i = 0; i < size; i++)
	{
		int u = delArray[i];
		for (int j = cd[u]; j < cd[u] + deg[u]; j++)
		{
			int v = adj[j];
			if (delMark[v] == true)
			{
				interEdges++;
				continue;
			}
			for (int k = cd[v]; k < cd[v] + deg[v]; k++)
			{
				int w = adj[k];
				if (w == u)
				{
					adj[k] = adj[cd[v] + deg[v] - 1];
					adj[cd[v] + deg[v] - 1] = w;
					deg[v]--;
					break;
				}
			}
		}
		deg[u] = 0;
	}
	delete[] delMark;
	return interEdges;
}

void Graph::deleteEdge(int s, int t)
{
	for (int j = cd[s]; j < cd[s] + deg[s]; j++)
	{
		int nbr = adj[j];

		if (nbr == t)
		{
			swap(adj[j], adj[cd[s] + deg[s] - 1]);
			deg[s]--;
			break;
		}
	}

	for (int j = cd[t]; j < cd[t] + deg[t]; j++)
	{
		int nbr = adj[j];

		if (nbr == s)
		{
			swap(adj[j], adj[cd[t] + deg[t] - 1]);
			deg[t]--;
			break;
		}
	}
}

void Graph::insertEdge(int s, int t)
{
	adj[cd[s] + deg[s]++] = t;
	adj[cd[t] + deg[t]++] = s;
}