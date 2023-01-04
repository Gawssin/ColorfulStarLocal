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
#include <cstdarg>
#include <unordered_map>
#include <omp.h>
#include "../header/heapLLU.h"
#include "../header/Graph.hpp"
#include "../header/tool.hpp"
#include "../header/ColorfulStarCore.hpp"

using namespace std;

int main(int argc, char **argv)
{


    auto t0 = getTime();
    Graph g;
    int h = atoi(argv[1]);
    cout << "Reading edgelist from file: " << argv[2] << endl;
    g.readedgelist(argv[2]);
    cout << "Reading edgelist finished!" << endl;
    g.mkGraph();
    cout << "mkGraph finished!" << endl;
    
    omp_set_num_threads(atoi(argv[3]));
    
    auto t1 = getTime();

    int *color = new int[g.n];
    int colorNum = g.color(color);
    
    printf("Total number of colors: %d\n", colorNum);


    __int128 **dp = new __int128 *[g.n];
    int **CC = new int *[g.n];
    __int128 *ColofulStarCoreNum = new __int128[g.n];

    initColStarDegree(g, dp, h, colorNum, color, CC);
    ColorfulStarHIndexSync(g, dp, h, color, CC, ColofulStarCoreNum, colorNum);


    auto t2 = getTime();
    printf("- Overall time = %lf\n", ((double)timeGap(t1, t2)) / 1e6);

    return 0;
}
