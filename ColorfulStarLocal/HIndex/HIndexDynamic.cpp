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
#include "../header/heapLLU.h"
#include "../header/Graph.hpp"
#include "../header/tool.hpp"
#include "../header/ColorfulStarCoreDynamic.hpp"

using namespace std;

struct dynamicEdge
{
    int s, t;
    bool status; // true for insertion, false for deletion
};

// boundBFS(g, qNbr, oldHindexS, newHindexS, numNbrBwn, vis);
void boundBFS(Graph &g, queue<int> &qNbr, __int128 oldHindexS, __int128 newHindexS, __int128 *ColofulStarCoreNum,
              int &numNbrBwn, int *affectedNodes, bool *vis)
{
    while (qNbr.empty() == false)
    {
        int curNode = qNbr.front();
        qNbr.pop();
        // affectedNodes[numNbrBwn++] = curNode;

        for (int j = g.cd[curNode]; j < g.cd[curNode] + g.deg[curNode]; j++)
        {
            int nbr = g.adj[j];

            if (vis[nbr] == false && ColofulStarCoreNum[nbr] <= oldHindexS && ColofulStarCoreNum[nbr] > newHindexS)
            {
                vis[nbr] = true;
                // numNbrBwn++;
                affectedNodes[numNbrBwn++] = nbr;
                qNbr.push(nbr);
            }
        }
    }
}

int main(int argc, char **argv)
{
    // cout << (1<<31)-1 << endl;
    // __int128 aaa = (__int128(1) << 127) - 1;
    // cout << _int128_to_str(aaa) << endl;
    // cout << _int128_to_str(numeric_limits<__int128>::max()) << endl;
    char *argv1, *argv2, *argv3;
    argv1 = argv[1], argv2 = argv[2];

    int numOfEdges = atoi(argv[3]);

    // int basicFlag = 0, colorAlgo = atoi(argv[3]);
    // if (argc >= 5 && strcmp(argv[4], "basic") == 0)
    // {
    //     basicFlag = 1;
    // }

    // cout << "Version: " << (basicFlag ? "Basic (HStarDP)" : "Advanced (HStarCD)") << endl << endl;

    auto t0 = getTime();

    Graph g;
    int h = atoi(argv1);
    cout << "Reading edgelist from file: " << argv2 << endl;
    g.readedgelist(argv2);
    cout << "Reading edgelist finished!" << endl;
    g.mkGraph();
    cout << "mkGraph finished!" << endl;
    auto t1 = getTime();

    int *edgeId = new int[g.e];
    for (int i = 0; i < g.e; i++)
        edgeId[i] = i;

    unsigned seed = system_clock::now().time_since_epoch().count();
    shuffle(edgeId, edgeId + g.e, default_random_engine(seed));

    struct dynamicEdge *dEdges = new struct dynamicEdge[numOfEdges];
    for (int i = 0; i < numOfEdges; i++)
    {
        dEdges[i].s = g.edges[edgeId[i]].s;
        dEdges[i].t = g.edges[edgeId[i]].t;
        dEdges[i].status = false;
        // cout << dEdges[i].s << " " << dEdges[i].t << endl;
    }

    // cout << g.isEdge(3,10) << endl;;
    //-----test
    //  numOfEdges = 1;
    //  dEdges[0].s = 0;
    //  dEdges[0].t = 1;
    //  dEdges[0].status = false;

    //-----test

    int *color = new int[g.n];
    int colorNum = g.color(color);
    printf("Total number of colors: %d\n", colorNum);

    // for (int i = 0; i < g.n; i++)
    // {
    //     printf("%d:(deg:%d) -> %d\n", i, g.deg[i], color[i]);
    // }

    __int128 **dp = new __int128 *[g.n];
    int **CC = new int *[g.n];
    __int128 *ColofulStarCoreNum = new __int128[g.n];

    //----------------------------batch deletion

    // for (int i = 0; i < numOfEdges; i++)
    // {
    //     g.deleteEdge(dEdges[i].s, dEdges[i].t);
    // }

    // initColStarDegree(g, dp, h, colorNum, color, CC);
    // ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);

    // if (argc >= 6)
    // {
    //     ofstream outfile(argv[5], ios::out);

    //     for (int i = 0; i < g.n; i++)
    //     {
    //         outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    //     }

    //     outfile.close();
    // }

    // for (int i = 0; i < numOfEdges; i++)
    // {
    //     g.insertEdge(dEdges[i].s, dEdges[i].t);
    // }

    //----------------------------batch deletion

    initColStarDegree(g, dp, h, colorNum, color, CC);

    ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);

    // for (int i = 0; i < g.n; i++)
    // {
    //     cout << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    // }

    // for (int i = 0; i < g.n; i++)
    // if(ColofulStarCoreNum[i] == 234136)
    // {
    //     cout << i << endl;
    //     return 0;
    // }

    // for (int i = 0; i < g.n; i++)
    // {
    //     if(ColofulStarCoreNum[i] < 6220 && ColofulStarCoreNum[i] > 6200)
    //         cout << i << " ++++++ " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    // }

    int ReComp = 0;

    if (argc >= 5)
    {
        ReComp = atoi(argv[4]);
    }

    __int128 *CoreNumBak = new __int128[g.n];
    for (int i = 0; i < g.n; i++)
    {
        CoreNumBak[i] = ColofulStarCoreNum[i];
    }

    auto t2 = getTime();
    if (ReComp == 1)
    {
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

        ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);

        // if (argc >= 5)
        // {
        //     ofstream outfile(argv[4], ios::out);

        //     for (int i = 0; i < g.n; i++)
        //     {
        //         outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
        //     }

        //     outfile.close();
        // }

        // return 0;

        for (int i = 0; i < numOfEdges; i++)
        {
            g.insertEdge(dEdges[i].s, dEdges[i].t);
        }
    }

    auto t3 = getTime();

    swap(CoreNumBak, ColofulStarCoreNum);

    auto t5 = getTime();

    cout << "insertion Total time (ms):\t" << timeGap(t3, t5) / 1000.0 << endl;

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
            // cout << "numNbrBwn: " << numNbrBwn << "\t" << cnt << endl;
        }
    }

    cout << "tolNumEstimated:\t" << tolNumEstimated << endl;
    cout << "tolNumChanged:\t" << tolNumChanged << endl;
    cout << "hit rate:\t" << 1.0 * tolNumChanged / tolNumEstimated * 100 << "%" << endl;

    cout << "ReComp Total time (ms):\t" << timeGap(t2, t3) / 1000.0 << endl;

    auto t4 = getTime();

    cout << "Time per edge (ms):\t" << timeGap(t5, t4) / 1000.0 / numOfEdges << endl;
    cout << "Total time (ms):\t" << timeGap(t5, t4) / 1000.0 << endl;

    // if (argc >= 5)
    // {
    //     ofstream outfile(argv[4], ios::out);

    //     for (int i = 0; i < g.n; i++)
    //     {
    //         outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    //     }

    //     outfile.close();
    // }

    return 0;

    //-----------------------------------
    // __int128 *coreTMP = new __int128[g.n];

    // memcpy(coreTMP, ColofulStarCoreNum, sizeof(__int128) * g.n);

    // int dId = atoi(argv[3]); // 10 --- 654

    // cout << dId << " ------ " << _int128_to_str(ColofulStarCoreNum[dId]) << endl;
    // for (int j = g.cd[dId]; j < g.cd[dId] + g.deg[dId]; j++)
    // {
    //     int nbr = g.adj[j];

    //     cout << nbr << " -> " << _int128_to_str(ColofulStarCoreNum[nbr]) << endl;
    // }

    // int tId = atoi(argv[4]);

    // __int128 oldHindex = ColofulStarCoreNum[tId];
    // __int128 newHindex = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, tId);

    // queue<int> qNbr;
    // qNbr.push(tId);

    // bool *vis = new bool[g.n]();
    // int numNbrBwn = 1;
    // vis[tId] = true;

    // while (qNbr.empty() == false)
    // {
    //     int curNode = qNbr.front();
    //     qNbr.pop();

    //     for (int j = g.cd[curNode]; j < g.cd[curNode] + g.deg[curNode]; j++)
    //     {
    //         int nbr = g.adj[j];

    //         if (vis[nbr] == false && ColofulStarCoreNum[nbr] <= oldHindex && ColofulStarCoreNum[nbr] > newHindex)
    //         {
    //             vis[nbr] = true;
    //             numNbrBwn++;
    //             qNbr.push(nbr);
    //         }
    //     }
    // }

    // cout << "number of nbrs between oldHindex and newHindex: " << numNbrBwn << endl;

    // dp[tId][h - 1] = oldHindex;
    // ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);

    // int innerNum = 0;
    // for (int i = 0; i < g.n; i++)
    // {
    //     if (ColofulStarCoreNum[i] != coreTMP[i])
    //     {

    //         cout << i << " : " << _int128_to_str(coreTMP[i]) << " --- " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    //     }

    //     if (ColofulStarCoreNum[i] <= oldHindex && ColofulStarCoreNum[i] >= newHindex)
    //         innerNum++;
    // }

    // cout << "innerNum: " << innerNum << endl;

    // // for (int i = 0; i < g.n; i++)
    // // {
    // //     printf("%d -> %s\n", i, _int128_to_str(dp[i][h - 1]));
    // // }

    // // ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, basicFlag);

    // auto t2 = getTime();
    // printf("- Overall time = %lfs\n", ((double)timeGap(t1, t2)) / 1e6);

    //-----------------------------------

    // // --------------------------------------
    // if (argc >= 3)
    // {
    //     ofstream outfile(argv[3], ios::out);

    //     for (int i = 0; i < g.n; i++)
    //     {
    //         outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
    //     }

    //     outfile.close();
    // }
    // // --------------------------------------

    // study nodes' distribution

    // if (argc >= 3)
    // {
    //     unordered_map<__int128, int> umap;
    //     for (int i = 0; i < g.n; i++)
    //     {
    //         umap[dp[i][h - 1]]++;
    //     }

    //     ofstream outfile(argv[3], ios::out);

    //     for (auto x : umap)
    //         outfile << _int128_to_str(x.first) << " " << x.second << endl;

    //     outfile.close();
    // }

    return 0;
}
