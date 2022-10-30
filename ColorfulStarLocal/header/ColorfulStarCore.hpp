
void initColStarDegree(Graph& g, __int128** dp, int h, int colorNum, int* color, int** CC, int basicFlag = 0)
{
    __int128* NotColor0 = new __int128[h]();
    __int128* MustColor0 = new __int128[h]();
    for (int i = 0; i < g.n; i++)
    {
        dp[i] = new __int128[h];
        int colorNum_i = 0;
        //int* C = new int[colorNum]();
        CC[i] = new int[colorNum]();
        for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
        {
            int nbr = g.adj[j];
            CC[i][color[nbr]]++;
            colorNum_i = max(colorNum_i, color[nbr]);
        }

        if (basicFlag == 0)
        {
            NotColor0[0] = 1;
            for (int c = 1; c <= colorNum_i; c++)
            {
                for (int j = h - 1; j > 0; j--)
                    NotColor0[j] = NotColor0[j - 1] * CC[i][c] + NotColor0[j];
            }

            for (int j = 1; j < h; j++)
            {
                MustColor0[j] = NotColor0[j - 1] * CC[i][0];
                dp[i][j] = MustColor0[j] + NotColor0[j];
                NotColor0[j - 1] = 0;
            }
            NotColor0[h - 1] = 0;
        }
        else
        {
            dp[i][0] = 1;
            for (int c = 0; c <= colorNum_i; c++)
            {
                for (int j = h - 1; j > 0; j--)
                    dp[i][j] = dp[i][j] + dp[i][j - 1] * CC[i][c];
            }
        }
    }
    delete[] NotColor0;
    delete[] MustColor0;
}


void ColorfulStarCoreDecomp(Graph& g, __int128** dp, int h, int* color, int** CC, __int128* ColofulStarCoreNum, int colorNum, int basicFlag, __int128* maxCore = 0, int* maxCoreNum = 0)
{
    __int128* tmpDP = new __int128[g.n];
    for (int i = 0; i < g.n; i++) tmpDP[i] = dp[i][h - 1];
    bheapLLU<__int128>* heap = mkheapLLU<__int128>(g.n, tmpDP);

    __int128 maxStarDegree = -1;
    for (int i = 0; i < g.n; i++)
    {
        maxStarDegree = max(maxStarDegree, dp[i][h - 1]);
    }
    printf("The maximum colorful h-star degree: %s\n", _int128_to_str(maxStarDegree));



    int leftN = g.n, leftM = g.e;

    __int128* NotColor = new __int128[h]();
    __int128* MustColor = new __int128[h]();

    __int128 starCoreNum = 0;
    int times = 0, maxN = 0, maxM = 0;
    keyvalueLLU<__int128> kv;

    printf("Times\t\tLeft_Nodes\tMaxCore(N)\tMaxCore(M)\tMaxCore(Density)\tMaxCore(CoreNumber)\n");

    while (leftN > 0)
    {
        times++;
        kv = popminLLU<__int128>(heap);

        if (kv.value > starCoreNum)
        {
            starCoreNum = kv.value;
            maxN = leftN;
            maxM = leftM;
        }

        //if (times % 100000 == 0)
            //printf("times = %d\tleft nodes = %-10d\ttolMax = %lf\tmaxN = %d\tmaxM = %d\tdensity = %lf\n", times, leftN, starCoreNum, maxN, maxM, maxN?(1.0 * maxM / maxN) : 0.0);

        if (times % 100000 == 0)
            printf("%d\t\t%d\t\t%d\t\t%d\t\t%lf\t\t%s\n", times, leftN, maxN, maxM, maxN ? (1.0 * maxM / maxN) : 0.0, _int128_to_str(starCoreNum));

        leftN--;
        int i = kv.key;
        leftM -= g.deg[i];
        ColofulStarCoreNum[i] = starCoreNum;

        for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
        {
            int nbr = g.adj[j];

            for (int p = g.cd[nbr]; p < g.cd[nbr] + g.deg[nbr]; p++)
            {
                int hnbr = g.adj[p];
                if (hnbr == i)
                {
                    swap(g.adj[p], g.adj[g.cd[nbr] + g.deg[nbr] - 1]);
                    g.deg[nbr]--;
                    break;
                }
            }

            CC[nbr][color[i]]--;

            NotColor[0] = 1;
            if (basicFlag == 0)
            {
                for (int p = 1; p < h; p++)
                {
                    MustColor[p] = NotColor[p - 1] * (CC[nbr][color[i]] + 1);
                    NotColor[p] = dp[nbr][p] - MustColor[p];
                }
                for (int p = 1; p < h; p++)
                {
                    MustColor[p] = NotColor[p - 1] * CC[nbr][color[i]];
                    dp[nbr][p] = NotColor[p] + MustColor[p];
                }
            }
            else
            {
                fill(dp[nbr] + 1, dp[nbr] + h, 0);
                for (int c = 0; c < colorNum; c++)
                {
                    for (int j = h - 1; j > 0; j--)
                        dp[nbr][j] = dp[nbr][j] + dp[nbr][j - 1] * CC[nbr][c];
                }
            }

            updateLLU(heap, nbr, dp[nbr][h - 1]);
        }
        g.deg[i] = 0;
    }

    if (maxCore != 0)*maxCore = starCoreNum;
    if (maxCoreNum != 0)*maxCoreNum = maxN;

    printf("End:\n%d\t\t%d\t\t%d\t\t%d\t\t%lf\t\t%s\n", times, leftN, maxN, maxM, maxN ? (1.0 * maxM / maxN) : 0.0, _int128_to_str(starCoreNum));

    delete[] NotColor;
    delete[] MustColor;
}

bool ubCMP(const pair<int, __int128>& a, const pair<int, __int128>& b)
{
    if (a.second == b.second)
        return a.first > b.first;
    return a.second > b.second;
}

void insertSort(pair<int, __int128>* idUB, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            int p = j + 1;
            if ((idUB[j].second == idUB[p].second && idUB[j].first < idUB[p].first) ||
                idUB[j].second < idUB[p].second)
            {
                swap(idUB[j].first, idUB[p].first);
                swap(idUB[j].second, idUB[p].second);
            }
            else
            {
                break;
            }
        }
    }
}

__int128 compHIndex(Graph& g, int u, __int128** dp, __int128* uBound, int h, __int128* nbrUB, int* color, int** CC, int colorNum,
    pair<int, __int128>* idUB, __int128* NotColor, __int128* MustColor, int* lastPos, int* lastId, __int128* boundCore, bool* prune)
{
    // pair<int, __int128>* idUB = new pair<int, __int128>[g.deg[u]];

    // __int128* NotColor = new __int128[h]();
    // __int128* MustColor = new __int128[h]();

    for (int i = g.cd[u]; i < g.cd[u] + g.deg[u]; i++)
    {
        int nbrID = g.adj[i];
        idUB[i - g.cd[u]].first = nbrID;
        idUB[i - g.cd[u]].second = uBound[nbrID];
    }

    if (lastId[u] == -1)
        sort(idUB, idUB + g.deg[u], ubCMP);
    else insertSort(idUB, g.deg[u]);

    for (int i = g.cd[u]; i < g.cd[u] + g.deg[u]; i++)
    {
        g.adj[i] = idUB[i - g.cd[u]].first;
    }

    // if (idUB[lastPos[u]].first == lastId[u] && idUB[lastPos[u]].second == boundCore[u])
    // {
    //     //printf("u = %d\n", u);
    //     *prune = true;
    //     return dp[u][h - 1];
    // }




    // for (int i = g.cd[u]; i < g.cd[u] + g.deg[u]; i++)
    // {
    //     int nbrID = g.adj[i];
    //     nbrUB[i - g.cd[u]] = uBound[nbrID];
    // }





    //sort(nbrUB, nbrUB + g.deg[u], ubCMP);

    // if (u == 1582)
    // {
    //     for (int i = 0; i < g.deg[u]; i++)
    //     {
    //         printf("u=1    %d -> %s\n", idUB[i].first, _int128_to_str(idUB[i].second));
    //     }

    // }

    // if (u == 313167)
    // {
    //     //printf("222 \n");
    //     for (int i = 0; i < g.deg[u]; i++)
    //     {
    //         printf("u=203032  color=%d,  %d -> %s\n", color[idUB[i].first], idUB[i].first, _int128_to_str(idUB[i].second));
    //     }

    // }





    // if (idUB[h - 2].second == 0)
    // {
    //     dp[u][h - 1] = 0;
    //     return 0;
    // }



    //__int128 lastUB = dp[u][h - 1];

    //fill(CC[u], CC[u] + colorNum, 0);
    memset(CC[u], 0, sizeof(int) * colorNum);
    int colorNum_u = 0;

    for (int j = 0; j < h - 2; j++)
    {
        int nbr = idUB[j].first;
        CC[u][color[nbr]]++;
        colorNum_u = max(colorNum_u, color[nbr]);
    }




    dp[u][0] = 1;
    //fill(dp[u] + 1, dp[u] + h, 0);
    memset(dp[u] + 1, 0, sizeof(__int128) * (h - 1));
    for (int c = 0; c <= colorNum_u; c++)
    {
        if (CC[u][c] > 0)
        {
            for (int j = h - 1; j > 0; j--)
                dp[u][j] = dp[u][j] + dp[u][j - 1] * CC[u][c];
        }
    }




    for (int L = h - 1; L <= g.deg[u]; L++)
    {
        int nbr = idUB[L - 1].first;
        CC[u][color[nbr]]++;
        colorNum_u = max(colorNum_u, color[nbr]);



        // for (int j = 0; j < L; j++)
        // {
        //     int nbr = idUB[j].first;
        //     CC[u][color[nbr]]++;
        //     colorNum_u = max(colorNum_u, color[nbr]);
        // }



        //===============================================

        // dp[u][0] = 1;
        // //fill(dp[u] + 1, dp[u] + h, 0);
        // memset(dp[u] + 1, 0, sizeof(__int128) * (h - 1));
        // for (int c = 0; c <= colorNum_u; c++)
        // {
        //     if (CC[u][c] > 0)
        //     {
        //         for (int j = h - 1; j > 0; j--)
        //             dp[u][j] = dp[u][j] + dp[u][j - 1] * CC[u][c];
        //     }
        // }

        //===============================================



        //-----------------------------------------------

        NotColor[0] = 1;
        //dp[u][0] = 1;
        //memset(dp[u] + 1, 0, sizeof(__int128) * (h - 1));

        for (int p = 1; p < h; p++)
        {
            MustColor[p] = NotColor[p - 1] * (CC[u][color[nbr]] - 1);
            NotColor[p] = dp[u][p] - MustColor[p];
        }
        for (int p = 1; p < h; p++)
        {
            MustColor[p] = NotColor[p - 1] * CC[u][color[nbr]];
            dp[u][p] = NotColor[p] + MustColor[p];
        }


        //-----------------------------------------------






        // if (u == 1)
        // {

        //     printf("u=1 dp    %s\n", _int128_to_str(dp[u][h - 1]));


        // }


        if (L == g.deg[u] || !(dp[u][h - 1] < idUB[L - 1].second && dp[u][h - 1] < idUB[L].second))
        {
            dp[u][h - 1] = min(dp[u][h - 1], idUB[L - 1].second);
            lastPos[u] = L - 1;
            boundCore[u] = idUB[L - 1].second;
            lastId[u] = idUB[L - 1].first;
            break;
        }


        //====================================

        // if (dp[u][h - 1] >= idUB[L - 1].second)  //dp[u][h - 1] = 1, idUB[L - 1].second = 0
        // {
        //     dp[u][h - 1] = idUB[L - 1].second;
        //     break;
        // }

        // if ((L < g.deg[u] && dp[u][h - 1] >= idUB[L].second) || L == g.deg[u]) break;

        //====================================







        // printf("Total number of colors\n");

        // if (L < g.deg[u])
        // {
        //     if (lastUB >= idUB[L].second)
        //         break;
        // }
        // else if (L == g.deg[u])
        //     break;

        //}

    }
    // if (u == 818)
    // {
    //     printf("333 \n");
    // }

    // delete[] idUB;
    // delete[] NotColor;
    // delete[] MustColor;

    return dp[u][h - 1];    //deg < L-1
}




void ColorfulStarHIndex(Graph& g, __int128** dp, int h, int* color, int** CC, __int128* ColofulStarCoreNum, int colorNum)
{
    bool updateFlag = true;
    bool* flag = new bool[g.n], * lastFlag = new bool[g.n];

    __int128* uBound = new __int128[g.n];
    __int128* nbrUB = new __int128[g.maxDeg];
    for (int i = 0; i < g.n; i++) uBound[i] = dp[i][h - 1], flag[i] = true, lastFlag[i] = true;

    int times = 0;



    pair<int, __int128>* idUB = new pair<int, __int128>[g.maxDeg];

    __int128* NotColor = new __int128[h]();
    __int128* MustColor = new __int128[h]();

    int* lastPos = new int[g.n]();
    int* lastId = new int[g.n]();

    fill(lastId, lastId + g.n, -1);
    __int128* boundCore = new __int128[g.n]();

    int invokeNum = 0, remainNum = 0;


    while (updateFlag)
    {
        printf("updating %d\n", times++);
        updateFlag = false;

        int unchangeNum = 0;

        invokeNum = 0, remainNum = 0;

        for (int i = 0; i < g.n; i++)
        {

            // if (i == 313167)
            // {
            //     printf("ub: %s\n", _int128_to_str(uBound[i]));
            //     for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
            //     {
            //         int nbr = g.adj[j];
            //         printf("u=313167  color=%d,  %d -> %s\n", color[nbr], nbr, _int128_to_str(uBound[nbr]));
            //     }
            // }

            if (uBound[i] == 0)
            {
                flag[i] = false;
                continue;
            }
            bool updatedi = false;
            for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
            {
                int nbrID = g.adj[j];

                if (lastId[i] == -1)
                {
                    updatedi = true;
                    break;
                }

                if (j - g.cd[i] > lastPos[i])
                    break;

                //if (lastFlag[nbrID] == true)
                if (lastFlag[nbrID] == true && uBound[nbrID] < uBound[i])
                {
                    updatedi = true;
                    break;
                }
            }

            if (updatedi == false) continue;

            //printf("--------------          %d -> %s\n", i, _int128_to_str(uBound[i]));



            __int128 ubBak = uBound[i];
            bool prune = false;

            uBound[i] = compHIndex(g, i, dp, uBound, h, nbrUB, color, CC, colorNum, idUB, NotColor, MustColor, lastPos, lastId, boundCore, &prune);

            invokeNum++;


            //printf("-------------- after    %d -> %s\n", i, _int128_to_str(uBound[i]));

            // if (ubBak < uBound[i])
            // {
            //     printf("larger:    %d -> %s --- %s\n", i, _int128_to_str(ubBak), _int128_to_str(uBound[i]));

            //     for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
            //     {
            //         int nbrID = g.adj[j];

            //         printf("%d -> %d\n", nbrID, color[nbrID]);

            //     }
            //     exit(0);
            // }
            
            // if(i == 85371)
            // {
            //     cout << _int128_to_str(ubBak) << " ooooooooooooooooo " << _int128_to_str(uBound[i]) << endl;
            // }


            if (ubBak != uBound[i])
            {
                updateFlag = true;
                flag[i] = true;
            }
            else
            {
                remainNum++;
                flag[i] = false;
                if (prune == true)
                    unchangeNum++;
            }

        }

        swap(flag, lastFlag);

        printf("unchangeNum: %d\n", unchangeNum);

        printf("invoke: %d\tremain: %d\n", invokeNum, remainNum);






        // for (int i = 0; i < g.n; i++)
        // {
        //     printf("--------------    %d -> %s\n", i, _int128_to_str(dp[i][h - 1]));
        // }
        //break;
        //-------------------------------------
    }


    for (int i = 0; i < g.n; i++)
    {
        ColofulStarCoreNum[i] = uBound[i];
        //printf("uBound: %d -> %s\n", i, _int128_to_str(uBound[i]));
    }

}


__int128 ColorfulStarHIndexSingel(Graph& g, __int128** dp, int h, int* color, int** CC, __int128* ColofulStarCoreNum, int colorNum, int nodeId)
{
    bool updateFlag = true;
    bool* flag = new bool[g.n], * lastFlag = new bool[g.n];

    __int128* uBound = new __int128[g.n];
    __int128* nbrUB = new __int128[g.maxDeg];
    for (int i = 0; i < g.n; i++) uBound[i] = dp[i][h - 1], flag[i] = true, lastFlag[i] = true;




    pair<int, __int128>* idUB = new pair<int, __int128>[g.maxDeg];

    __int128* NotColor = new __int128[h]();
    __int128* MustColor = new __int128[h]();

    int* lastPos = new int[g.n]();
    int* lastId = new int[g.n]();

    fill(lastId, lastId + g.n, -1);
    __int128* boundCore = new __int128[g.n]();

    
    bool prune = false;

    __int128 imHindex = compHIndex(g, nodeId, dp, uBound, h, nbrUB, color, CC, colorNum, idUB, NotColor, MustColor, lastPos, lastId, boundCore, &prune);

    cout << nodeId << " -- ---- now ---- -- " << _int128_to_str(imHindex) << endl;
    return imHindex;
}


