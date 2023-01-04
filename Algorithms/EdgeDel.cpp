#include <stdlib.h>
#include <cstdio>
#include <stdbool.h>
#include <string.h>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <queue>
#include <cstdarg>
#include <random>
#include <unordered_map>
#include <omp.h>
#include "../header/heapLLU.h"
#include "../header/Graph.hpp"
#include "../header/tool.hpp"
#include "../header/ColorfulStarCore.hpp"

using namespace std;

struct dynamicEdge
{
    int s, t;
    bool status; // true for insertion, false for deletion
};

void boundBFS(Graph &g, queue<int> &qNbr, __int128 oldHindexS, __int128 newHindexS, __int128 *ColofulStarCoreNum,
              int &numNbrBwn, int *affectedNodes, bool *vis)
{
    while (qNbr.empty() == false)
    {
        int curNode = qNbr.front();
        qNbr.pop();

        for (int j = g.cd[curNode]; j < g.cd[curNode] + g.deg[curNode]; j++)
        {
            int nbr = g.adj[j];

            if (vis[nbr] == false && ColofulStarCoreNum[nbr] <= oldHindexS && ColofulStarCoreNum[nbr] > newHindexS)
            {
                vis[nbr] = true;
                affectedNodes[numNbrBwn++] = nbr;
                qNbr.push(nbr);
            }
        }
    }
}

int main(int argc, char **argv)
{
    auto t0 = getTime();
    omp_set_num_threads(atoi(argv[5]));

    Graph g;
    int h = atoi(argv[1]);
    cout << "Reading edgelist from file: " << argv[2] << endl;
    g.readedgelist(argv[2]);
    cout << "Reading edgelist finished!" << endl;
    g.mkGraph();
    cout << "mkGraph finished!" << endl;
    auto t1 = getTime();

    int numOfEdges = atoi(argv[3]);
    int *edgeId = new int[g.e];
    for (int i = 0; i < g.e; i++)
        edgeId[i] = i;

    unsigned seed = system_clock::now().time_since_epoch().count();
    shuffle(edgeId, edgeId + g.e, default_random_engine(seed));

    // randomly select edges with the number equal to "numOfEdges"
    struct dynamicEdge *dEdges = new struct dynamicEdge[numOfEdges];
    for (int i = 0; i < numOfEdges; i++)
    {
        dEdges[i].s = g.edges[edgeId[i]].s;
        dEdges[i].t = g.edges[edgeId[i]].t;
        dEdges[i].status = false;
    }

    int *color = new int[g.n];
    int colorNum = g.color(color);
    printf("Total number of colors: %d\n", colorNum);

    __int128 **dp = new __int128 *[g.n];
    int **CC = new int *[g.n];
    __int128 *ColofulStarCoreNum = new __int128[g.n];

    int optimization = 3;
    initColStarDegree(g, dp, h, colorNum, color, CC);
    ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, optimization);

    int ReComp = 0;

    if (argc >= 5)
    {
        ReComp = atoi(argv[4]);
    }

    auto t2 = getTime();
    if (ReComp == 1)
    {
        __int128 *CoreNumBak = new __int128[g.n];
        for (int i = 0; i < g.n; i++)
        {
            CoreNumBak[i] = ColofulStarCoreNum[i];
        }
        for (int processNum = 0; processNum < numOfEdges; processNum++)
        {
            int sId = dEdges[processNum].s;
            int tId = dEdges[processNum].t;
            if (dEdges[processNum].status == false)
            {
                g.deleteEdge(sId, tId);
            }
        }

        for (int i = 0; i < g.n; i++)
            dp[i][h - 1] = ColofulStarCoreNum[i];

        ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, optimization);

        for (int i = 0; i < numOfEdges; i++)
        {
            g.insertEdge(dEdges[i].s, dEdges[i].t);
        }
        swap(CoreNumBak, ColofulStarCoreNum);
    }

    auto t4 = getTime();

    queue<int> qNbr;
    bool *vis = new bool[g.n]();
    int numNbrBwn;
    int *affectedNodes = new int[g.n];

    __int128 *coreBak = new __int128[g.n];

    long long tolNumEstimated = 0, tolNumChanged = 0;

    for (int processNum = 0; processNum < numOfEdges; processNum++)
    {
        int sId = dEdges[processNum].s;
        int tId = dEdges[processNum].t;

        if (dEdges[processNum].status == false)
        {
            g.deleteEdge(sId, tId);

            __int128 oldHindexS = ColofulStarCoreNum[sId];
            __int128 oldHindexT = ColofulStarCoreNum[tId];

            __int128 newHindexS, newHindexT;

            if (oldHindexS != oldHindexT) // if olds > news, news are the same?
            {
                if (oldHindexS > oldHindexT)
                {
                    swap(sId, tId);
                    swap(oldHindexS, oldHindexT);
                }
                newHindexS = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, sId);

                qNbr.push(sId);
                numNbrBwn = 0;
                memset(vis, 0, sizeof(bool) * g.n);
                vis[sId] = true;
                ColofulStarCoreNum[sId] = newHindexS;

                boundBFS(g, qNbr, oldHindexS, newHindexS, ColofulStarCoreNum, numNbrBwn, affectedNodes, vis);
            }
            else
            {
                newHindexS = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, sId);
                newHindexT = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, tId);

                if (newHindexS > newHindexT)
                {
                    swap(sId, tId);
                    swap(newHindexS, newHindexT);
                }

                qNbr.push(sId);
                qNbr.push(tId);
                numNbrBwn = 0;
                memset(vis, 0, sizeof(bool) * g.n);
                vis[sId] = true;
                vis[tId] = true;
                ColofulStarCoreNum[sId] = newHindexS;
                ColofulStarCoreNum[tId] = newHindexT;

                boundBFS(g, qNbr, oldHindexS, min(newHindexS, newHindexT), ColofulStarCoreNum, numNbrBwn, affectedNodes, vis);
                affectedNodes[numNbrBwn++] = tId;
            }

            for (int i = 0; i < numNbrBwn; i++)
            {
                coreBak[affectedNodes[i]] = ColofulStarCoreNum[affectedNodes[i]];
            }

            ColorfulStarHIndexNodes(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, numNbrBwn, affectedNodes);

            int cnt = 0;
            for (int i = 0; i < numNbrBwn; i++)
            {
                if (coreBak[affectedNodes[i]] != ColofulStarCoreNum[affectedNodes[i]])
                    cnt++;
            }

            tolNumEstimated += numNbrBwn;
            tolNumChanged += cnt;
        }
    }

    cout << "tolNumEstimated:\t" << tolNumEstimated << endl;
    cout << "tolNumChanged:\t" << tolNumChanged << endl;
    cout << "hit rate:\t" << 1.0 * tolNumChanged / tolNumEstimated * 100 << "%" << endl;

    cout << "ReComp Total time (ms):\t" << timeGap(t2, t4) / 1000.0 << endl;

    auto t5 = getTime();

    cout << "Time per edge (ms):\t" << timeGap(t4, t5) / 1000.0 / numOfEdges << endl;
    cout << "Total time (ms):\t" << timeGap(t4, t5) / 1000.0 << endl;

    // if (argc >= 6)
    // {
    //     ofstream outfile(argv[5], ios::out);

    //     for (int i = 0; i < g.n; i++)
    //     {
    //         outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    //     }

    //     outfile.close();
    // }

    return 0;
}
