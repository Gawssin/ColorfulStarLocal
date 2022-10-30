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

bool isClique(Graph &g, int numNbrBwn, int *affectedNodes)
{
    for (int i = 0; i < numNbrBwn; i++)
    {
        for (int j = i+1; j < numNbrBwn; j++)
        {
            if(g.isEdge(affectedNodes[i], affectedNodes[j]) == false)
                return false;
        }
    }
    return true;
}

// boundBFS(g, qNbr, oldHindexS, newHindexS, numNbrBwn, vis);
void boundBFS(Graph &g, queue<int> &qNbr, __int128 oldHindexS, __int128 newHindexS, __int128 *ColofulStarCoreNum,
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
                // numNbrBwn++;
                // affectedNodes[numNbrBwn++] = nbr;
                qNbr.push(nbr);
            }
        }
    }
}

int main(int argc, char **argv)
{
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
        dEdges[i].status = true;
        // cout << dEdges[i].s << " " << dEdges[i].t << endl;
    }

    // cout << g.isEdge(3,10) << endl;;
    //-----test

    //  numOfEdges = 1;
    //  dEdges[0].s = 454328  ;
    //  dEdges[0].t = 591165;
    //  dEdges[0].status = true;

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

    for (int i = 0; i < numOfEdges; i++)
    {
        g.deleteEdge(dEdges[i].s, dEdges[i].t);
    }

    initColStarDegree(g, dp, h, colorNum, color, CC);
    ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);
    
    __int128 maxcore = 0;
    for (int i = 0; i < g.n; i++)
    {
        maxcore = max(maxcore, ColofulStarCoreNum[i]);
    }
    cout << "maxcore: " << _int128_to_str(maxcore) << endl;
    

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

    // initColStarDegree(g, dp, h, colorNum, color, CC);

    // ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);

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
    int DynamicIns = 1;
    ReComp = atoi(argv[4]);

    if (argc >= 6)
    {
        DynamicIns = atoi(argv[5]);
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
                g.insertEdge(sId, tId);
            }
        }

        initColStarDegree(g, dp, h, colorNum, color, CC);
        ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);

        for (int i = 0; i < numOfEdges; i++)
        {
            g.deleteEdge(dEdges[i].s, dEdges[i].t);
        }
    }

    auto t3 = getTime();
    cout << "ReComp Time (ms):\t" << timeGap(t2, t3) / 1000.0 << endl;

    if (DynamicIns == 0)
        return 0;

    swap(CoreNumBak, ColofulStarCoreNum);

    queue<int> qNbr;
    bool *vis = new bool[g.n]();
    int numNbrBwn;
    int *affectedNodes = new int[g.n];

    __int128 *coreBak = new __int128[g.n];

    long long tolNumEstimated = 0, tolNumChanged = 0;

    auto t4 = getTime();
    
    // __int128 *starDegBak = new __int128[g.maxDeg];
    pair<int, __int128> *starDegBak = new pair<int, __int128>[g.maxDeg];
    for (int processNum = 0; processNum < numOfEdges; processNum++)
    {
        int sId = dEdges[processNum].s;
        int tId = dEdges[processNum].t;

        // cout << _int128_to_str(ColofulStarCoreNum[sId]) << " " << _int128_to_str(ColofulStarCoreNum[tId]) << endl;
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
                
                //[i - g.cd[sId]] = ColofulStarCoreNum[nbr];
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
            
            cout << "newHindexS:\t" <<  _int128_to_str(newHindexS) << endl;
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
            
            cout << "newHindexT:\t" <<  _int128_to_str(newHindexT) << endl;
            ColofulStarCoreNum[tId] = oldHindexT;
            
            

            //__int128 newHindexS = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, sId);

            // ColofulStarCoreNum[sId] = oldHindexS;

            __int128 UpBoundS, UpBoundT;

            // newHindexS = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, sId);

            ColofulStarCoreNum[tId] = (__int128(1) << 127) - 1;
            UpBoundS = ColorStarDegreeWithinCore(g, dp, h, colorNum, color, CC, ColofulStarCoreNum, sId);
            ColofulStarCoreNum[tId] = oldHindexT;

            ColofulStarCoreNum[sId] = (__int128(1) << 127) - 1;
            UpBoundT = ColorStarDegreeWithinCore(g, dp, h, colorNum, color, CC, ColofulStarCoreNum, tId);
            // newHindexT = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, tId);
            ColofulStarCoreNum[sId] = oldHindexS;

            __int128 upperBound = min(UpBoundS, UpBoundT);
            if (upperBound <= oldHindexS)
                continue;

            // if (oldHindexS != oldHindexT)
            // {
            //     __int128 newHindexS = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, sId);
            //     ColofulStarCoreNum[sId] = oldHindexS;

            //     if(newHindexS == oldHindexS)
            //     {
            //         cout << sId << "\t ------------------- \t" << tId << endl;
            //         //continue;

            //     }
            // }

            // if (upperBound < oldHindexT)

            qNbr.push(sId);
            numNbrBwn = 0;
            memset(vis, 0, sizeof(bool) * g.n);
            vis[sId] = true;
            // ColofulStarCoreNum[sId] = upperBound;

            boundBFS(g, qNbr, upperBound, oldHindexS, ColofulStarCoreNum, numNbrBwn, affectedNodes, vis);

            cout << "---> " << _int128_to_str(UpBoundS) << " " << _int128_to_str(UpBoundT) << " " 
                 << _int128_to_str(oldHindexS) << " " << numNbrBwn << endl;

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

            //static int* changed = new int[g.n];
            int cnt = 0;
            cout << sId << " " << tId << endl;
            for (int i = 0; i < numNbrBwn; i++)
            {
                if (coreBak[affectedNodes[i]] != ColofulStarCoreNum[affectedNodes[i]])
                {
                    if (coreBak[affectedNodes[i]] > ColofulStarCoreNum[affectedNodes[i]])
                    {
                        cout << "error ---------------" << endl;
                    }
                    
                    cnt++;
                    //changed[cnt-1] = affectedNodes[i];
                    maxCoreChanged = max(maxCoreChanged, ColofulStarCoreNum[affectedNodes[i]]);
                    //cout << _int128_to_str(coreBak[affectedNodes[i]]) << " " << _int128_to_str(ColofulStarCoreNum[affectedNodes[i]]) << endl;
                    //  if(coreBak[affectedNodes[i]] != oldHindexS)
                    //  {
                    //  cout << "------------------ " << _int128_to_str(ColofulStarCoreNum[affectedNodes[i]]) << endl;

                    // }
                }
            }
            
            // if(cnt > 0 && oldHindexS == oldHindexT && ColofulStarCoreNum[sId] == ColofulStarCoreNum[tId])
            // {
            //     cout << isClique(g, cnt, changed) << " -------------------------------------- Equal: " << _int128_to_str(oldHindexS) << "\t" << _int128_to_str(ColofulStarCoreNum[sId]) << endl;
                
            // }
            
            // tolNumEstimated += numNbrBwn;
            // tolNumChanged += cnt;

            // cout << _int128_to_str(oldHindexS) << "\t" << _int128_to_str(oldHindexT) << "\t" <<  _int128_to_str(upperBound) << "\t" << _int128_to_str(maxCoreChanged)
            //      << "\t" << g.deg[sId]<< endl;
           
            cout << _int128_to_str(oldHindexS) << "\t" << _int128_to_str(oldHindexT) << "\t numNbrBwn: " << numNbrBwn << "\t" << cnt << "\t" << _int128_to_str(maxCoreChanged) << endl;
        }
    }

    auto t5 = getTime();
    cout << "ReComp Time (ms):\t" << timeGap(t2, t3) / 1000.0 << endl;
    cout << "Time per edge (ms):\t" << timeGap(t4, t5) / 1000.0 / numOfEdges << endl;
    cout << "Total time (ms):\t" << timeGap(t4, t5) / 1000.0 << endl;

    // cout << "tolNumEstimated:\t" << tolNumEstimated << endl;
    // cout << "tolNumChanged:\t" << tolNumChanged << endl;
    // cout << "hit rate:\t" << 1.0 * tolNumChanged / tolNumEstimated * 100 << "%" << endl;

    if (argc >= 7)
    {
        ofstream outfile(argv[6], ios::out);

        for (int i = 0; i < g.n; i++)
        {
            outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
        }

        outfile.close();
    }

    return 0;
}
