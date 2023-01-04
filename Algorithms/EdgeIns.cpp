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
#include <omp.h>
#include <unordered_map>
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

bool isClique(Graph &g, int numNbrBwn, int *affectedNodes)
{
    for (int i = 0; i < numNbrBwn; i++)
    {
        for (int j = i + 1; j < numNbrBwn; j++)
        {
            if (g.isEdge(affectedNodes[i], affectedNodes[j]) == false)
                return false;
        }
    }
    return true;
}

void boundBFSIns(Graph &g, queue<int> &qNbr, __int128 oldHindexS, __int128 newHindexS, __int128 *ColofulStarCoreNum,
              int &numNbrBwn, int *affectedNodes, bool *vis)
{
    while (qNbr.empty() == false)
    {
        int curNode = qNbr.front();
        qNbr.pop();
        affectedNodes[numNbrBwn++] = curNode;

        for (int j = g.cd[curNode]; j < g.cd[curNode] + g.deg[curNode]; j++)
        {
            int nbr = g.adj[j];

            if (vis[nbr] == false && ColofulStarCoreNum[nbr] < oldHindexS && ColofulStarCoreNum[nbr] >= newHindexS)
            {
                vis[nbr] = true;
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
        dEdges[i].status = true;
    }

    int *color = new int[g.n];
    int colorNum = g.color(color);
    printf("Total number of colors: %d\n", colorNum);

    __int128 **dp = new __int128 *[g.n];
    int **CC = new int *[g.n];
    __int128 *ColofulStarCoreNum = new __int128[g.n];

    for (int i = 0; i < numOfEdges; i++)
    {
        g.deleteEdge(dEdges[i].s, dEdges[i].t);
    }
    int optimization = 3;

    initColStarDegree(g, dp, h, colorNum, color, CC);
    ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, optimization);

    __int128 maxcore = 0;
    for (int i = 0; i < g.n; i++)
    {
        maxcore = max(maxcore, ColofulStarCoreNum[i]);
    }
    cout << "maxcore: " << _int128_to_str(maxcore) << endl;

    int ReComp = 0;
    ReComp = atoi(argv[4]);

    __int128 *CoreNumBak = new __int128[g.n];

    auto t2 = getTime();
    if (ReComp == 1)
    {
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
                g.insertEdge(sId, tId);
            }
        }

        initColStarDegree(g, dp, h, colorNum, color, CC);
        ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, optimization);

        for (int i = 0; i < numOfEdges; i++)
        {
            g.deleteEdge(dEdges[i].s, dEdges[i].t);
        }
        swap(CoreNumBak, ColofulStarCoreNum);
    }

    auto t3 = getTime();

    queue<int> qNbr;
    bool *vis = new bool[g.n]();
    int numNbrBwn;
    int *affectedNodes = new int[g.n];

    __int128 *coreBak = new __int128[g.n];

    long long tolNumEstimated = 0, tolNumChanged = 0;

    auto t4 = getTime();

    pair<int, __int128> *starDegBak = new pair<int, __int128>[g.maxDeg];
    for (int processNum = 0; processNum < numOfEdges; processNum++)
    {
        int sId = dEdges[processNum].s;
        int tId = dEdges[processNum].t;

        __int128 maxCoreChanged = 0;
        if (dEdges[processNum].status == true)
        {
            g.insertEdge(sId, tId);
            __int128 oldHindexS = ColofulStarCoreNum[sId];
            __int128 oldHindexT = ColofulStarCoreNum[tId];

            if (oldHindexS > oldHindexT)
            {
                swap(sId, tId);
                swap(oldHindexS, oldHindexT);
            }
            __int128 lowerBound = oldHindexS;
            numNbrBwn = 0;

            for (int i = g.cd[sId]; i < g.cd[sId] + g.deg[sId]; i++)
            {
                int nbr = g.adj[i];
                starDegBak[i - g.cd[sId]].first = nbr;
                starDegBak[i - g.cd[sId]].second = ColofulStarCoreNum[nbr];

                if (ColofulStarCoreNum[nbr] >= lowerBound)
                {
                    affectedNodes[numNbrBwn++] = nbr;
                }
            }

            ColStarDegreeNodes(g, dp, h, colorNum, color, CC, ColofulStarCoreNum, numNbrBwn, affectedNodes, lowerBound);

            __int128 newHindexS = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, sId);

            for (int i = g.cd[sId]; i < g.cd[sId] + g.deg[sId]; i++)
            {
                g.adj[i] = starDegBak[i - g.cd[sId]].first;
                ColofulStarCoreNum[g.adj[i]] = starDegBak[i - g.cd[sId]].second;
            }

            // cout << "newHindexS:\t" <<  _int128_to_str(newHindexS) << endl;
            ColofulStarCoreNum[sId] = lowerBound;

            numNbrBwn = 0;
            for (int i = g.cd[tId]; i < g.cd[tId] + g.deg[tId]; i++)
            {
                int nbr = g.adj[i];
                starDegBak[i - g.cd[tId]].first = nbr;
                starDegBak[i - g.cd[tId]].second = ColofulStarCoreNum[nbr];
                if (ColofulStarCoreNum[nbr] >= lowerBound)
                {
                    affectedNodes[numNbrBwn++] = nbr;
                }
            }

            ColStarDegreeNodes(g, dp, h, colorNum, color, CC, ColofulStarCoreNum, numNbrBwn, affectedNodes, lowerBound);

            __int128 newHindexT = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, tId);

            for (int i = g.cd[tId]; i < g.cd[tId] + g.deg[tId]; i++)
            {
                g.adj[i] = starDegBak[i - g.cd[tId]].first;
                ColofulStarCoreNum[g.adj[i]] = starDegBak[i - g.cd[tId]].second;
            }

            // cout << "newHindexT:\t" <<  _int128_to_str(newHindexT) << endl;
            ColofulStarCoreNum[tId] = oldHindexT;

            __int128 UpBoundS, UpBoundT;

            ColofulStarCoreNum[tId] = (__int128(1) << 127) - 1;
            UpBoundS = ColorStarDegreeWithinCore(g, dp, h, colorNum, color, CC, ColofulStarCoreNum, sId);
            ColofulStarCoreNum[tId] = oldHindexT;

            ColofulStarCoreNum[sId] = (__int128(1) << 127) - 1;
            UpBoundT = ColorStarDegreeWithinCore(g, dp, h, colorNum, color, CC, ColofulStarCoreNum, tId);
            ColofulStarCoreNum[sId] = oldHindexS;

            __int128 upperBound = min(UpBoundS, UpBoundT);
            if (upperBound <= oldHindexS)
                continue;

            qNbr.push(sId);
            numNbrBwn = 0;
            memset(vis, 0, sizeof(bool) * g.n);
            vis[sId] = true;

            boundBFSIns(g, qNbr, upperBound, oldHindexS, ColofulStarCoreNum, numNbrBwn, affectedNodes, vis);

            // cout << "---> " << _int128_to_str(UpBoundS) << " " << _int128_to_str(UpBoundT) << " "
            //      << _int128_to_str(oldHindexS) << " " << numNbrBwn << endl;

            for (int i = 0; i < numNbrBwn; i++)
            {
                coreBak[affectedNodes[i]] = ColofulStarCoreNum[affectedNodes[i]];
            }

            ColStarDegreeNodes(g, dp, h, colorNum, color, CC, ColofulStarCoreNum, numNbrBwn, affectedNodes, oldHindexS);

            for (int i = 0; i < numNbrBwn; i++)
            {
                ColofulStarCoreNum[affectedNodes[i]] = min(upperBound, ColofulStarCoreNum[affectedNodes[i]]);
            }

            ColorfulStarHIndexNodes(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, numNbrBwn, affectedNodes);

            // int cnt = 0;
            // // cout << sId << " " << tId << endl;
            // for (int i = 0; i < numNbrBwn; i++)
            // {
            //     if (coreBak[affectedNodes[i]] != ColofulStarCoreNum[affectedNodes[i]])
            //     {
            //         if (coreBak[affectedNodes[i]] > ColofulStarCoreNum[affectedNodes[i]])
            //         {
            //             cout << "error." << endl;
            //         }
            //         cnt++;
            //         maxCoreChanged = max(maxCoreChanged, ColofulStarCoreNum[affectedNodes[i]]);
            //     }
            // }
            // cout << _int128_to_str(oldHindexS) << "\t" << _int128_to_str(oldHindexT) << "\t numNbrBwn: " << numNbrBwn << "\t" << cnt << "\t" << _int128_to_str(maxCoreChanged) << endl;
        }
    }

    auto t5 = getTime();
    cout << "ReComp Time (ms):\t" << timeGap(t2, t3) / 1000.0 << endl;
    cout << "Time per edge (ms):\t" << timeGap(t4, t5) / 1000.0 / numOfEdges << endl;
    cout << "Total time (ms):\t" << timeGap(t4, t5) / 1000.0 << endl;

    // if (argc >= 7)
    // {
    //     ofstream outfile(argv[6], ios::out);

    //     for (int i = 0; i < g.n; i++)
    //     {
    //         outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    //     }

    //     outfile.close();
    // }

    return 0;
}
