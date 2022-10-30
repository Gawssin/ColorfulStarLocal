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
#include <unordered_map>
#include "../header/heapLLU.h"
#include "../header/Graph.hpp"
#include "../header/tool.hpp"
#include "../header/ColorfulStarCore.hpp"

using namespace std;

int main(int argc, char **argv)
{
    char *argv1, *argv2, *argv3;
    argv1 = argv[1], argv2 = argv[2];
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

    initColStarDegree(g, dp, h, colorNum, color, CC);

    ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);
    
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

    __int128 *coreTMP = new __int128[g.n];

    memcpy(coreTMP, ColofulStarCoreNum, sizeof(__int128) * g.n);

    int dId = atoi(argv[3]); // 10 --- 654

    cout << dId << " ------ " << _int128_to_str(ColofulStarCoreNum[dId]) << endl;
    for (int j = g.cd[dId]; j < g.cd[dId] + g.deg[dId]; j++)
    {
        int nbr = g.adj[j];

        cout << nbr << " -> " << _int128_to_str(ColofulStarCoreNum[nbr]) << endl;
    }

    int tId = atoi(argv[4]);

    for (int j = g.cd[dId]; j < g.cd[dId] + g.deg[dId]; j++)
    {
        int nbr = g.adj[j];

        if (nbr == tId)
        {
            swap(g.adj[j], g.adj[g.cd[dId] + g.deg[dId] - 1]);
            g.deg[dId]--;
            break;
        }
    }

    for (int j = g.cd[tId]; j < g.cd[tId] + g.deg[tId]; j++)
    {
        int nbr = g.adj[j];

        if (nbr == dId)
        {
            swap(g.adj[j], g.adj[g.cd[tId] + g.deg[tId] - 1]);
            g.deg[tId]--;
            break;
        }
    }

    __int128 oldHindex = ColofulStarCoreNum[tId];
    __int128 newHindex = ColorfulStarHIndexSingel(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, tId);

    queue<int> qNbr;
    qNbr.push(tId);

    bool *vis = new bool[g.n]();
    int numNbrBwn = 1;
    vis[tId] = true;

    while (qNbr.empty() == false)
    {
        int curNode = qNbr.front();
        qNbr.pop();

        for (int j = g.cd[curNode]; j < g.cd[curNode] + g.deg[curNode]; j++)
        {
            int nbr = g.adj[j];

            if (vis[nbr] == false && ColofulStarCoreNum[nbr] <= oldHindex && ColofulStarCoreNum[nbr] > newHindex)
            {
                vis[nbr] = true;
                numNbrBwn++;
                qNbr.push(nbr);
            }
        }
    }
    
    cout << "number of nbrs between oldHindex and newHindex: " << numNbrBwn << endl;
    
    
    
    dp[tId][h-1] = oldHindex;
    ColorfulStarHIndex(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);

    int innerNum = 0;
    for (int i = 0; i < g.n; i++)
    {
        if (ColofulStarCoreNum[i] != coreTMP[i])
        {

            cout << i << " : " << _int128_to_str(coreTMP[i]) << " --- " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
        }

        if (ColofulStarCoreNum[i] <= oldHindex && ColofulStarCoreNum[i] >= newHindex)
            innerNum++;
    }

    cout << "innerNum: " << innerNum << endl;

    // for (int i = 0; i < g.n; i++)
    // {
    //     printf("%d -> %s\n", i, _int128_to_str(dp[i][h - 1]));
    // }

    // ColorfulStarCoreDecomp(g, dp, h, color, CC, ColofulStarCoreNum, colorNum, basicFlag);

    auto t2 = getTime();
    printf("- Overall time = %lfs\n", ((double)timeGap(t1, t2)) / 1e6);

    // // --------------------------------------
    if (argc >= 6)
    {
        ofstream outfile(argv[5], ios::out);

        for (int i = 0; i < g.n; i++)
        {
            outfile << i << " " << _int128_to_str(ColofulStarCoreNum[i]) << endl;
        }

        outfile.close();
    }
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
