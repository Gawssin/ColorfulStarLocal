
void initColStarDegree(Graph &g, __int128 **dp, int h, int colorNum, int *color, int **CC, int basicFlag = 1)
{
    __int128 *NotColor0 = new __int128[h]();
    __int128 *MustColor0 = new __int128[h]();

    for (int i = 0; i < g.n; i++)
    {
        dp[i] = new __int128[h];
        dp[i][0] = 1;
        int colorNum_i = 0;
        // int* C = new int[colorNum]();
        CC[i] = new int[colorNum]();
        for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
        {
            int nbr = g.adj[j];
            CC[i][color[nbr]]++;
            colorNum_i = max(colorNum_i, color[nbr]);
        }

        if (basicFlag == 0) // basicFlag = 1 if running parallelly
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

void ColorfulStarCoreDecomp(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 *ColofulStarCoreNum, int colorNum, int basicFlag, __int128 *maxCore = 0, int *maxCoreNum = 0)
{
    __int128 *tmpDP = new __int128[g.n];
    for (int i = 0; i < g.n; i++)
        tmpDP[i] = dp[i][h - 1];
    bheapLLU<__int128> *heap = mkheapLLU<__int128>(g.n, tmpDP);

    __int128 maxStarDegree = -1;
    for (int i = 0; i < g.n; i++)
    {
        maxStarDegree = max(maxStarDegree, dp[i][h - 1]);
    }
    printf("The maximum colorful h-star degree: %s\n", _int128_to_str(maxStarDegree));

    int leftN = g.n, leftM = g.e;

    __int128 *NotColor = new __int128[h]();
    __int128 *MustColor = new __int128[h]();

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

        // if (times % 100000 == 0)
        // printf("times = %d\tleft nodes = %-10d\ttolMax = %lf\tmaxN = %d\tmaxM = %d\tdensity = %lf\n", times, leftN, starCoreNum, maxN, maxM, maxN?(1.0 * maxM / maxN) : 0.0);

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

    if (maxCore != 0)
        *maxCore = starCoreNum;
    if (maxCoreNum != 0)
        *maxCoreNum = maxN;

    printf("End:\n%d\t\t%d\t\t%d\t\t%d\t\t%lf\t\t%s\n", times, leftN, maxN, maxM, maxN ? (1.0 * maxM / maxN) : 0.0, _int128_to_str(starCoreNum));

    delete[] NotColor;
    delete[] MustColor;
}

void ColorfulStarCore(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 LB, int &delNum, int &delEdges)
{
    printf("Lower Bound: %s\n", _int128_to_str(LB));

    __int128 *tmpDP = new __int128[g.n];
    for (int i = 0; i < g.n; i++)
        tmpDP[i] = dp[i][h - 1];

    bheapLLU<__int128> *heap = mkheapLLU<__int128>(g.n, tmpDP);

    __int128 maxStarDegree = -1;
    for (int i = 0; i < g.n; i++)
    {
        maxStarDegree = max(maxStarDegree, dp[i][h - 1]);
    }
    printf("maxStarDegree = %s\n", _int128_to_str(maxStarDegree));

    int leftN = g.n, leftM = g.e;

    __int128 *NotColor = new __int128[h]();
    __int128 *MustColor = new __int128[h]();

    __int128 starCoreNum = 0;
    int times = 0, maxN = 0, maxM = 0;
    keyvalueLLU<__int128> kv;

    while (leftN > 0)
    {
        times++;
        //__int128 Min = 1e300;
        kv = popminLLU<__int128>(heap);

        if (kv.value >= LB)
        {
            break;
        }
        delNum++;
        delEdges += g.deg[kv.key];

        if (kv.value > starCoreNum)
        {
            starCoreNum = kv.value;
            maxN = leftN;
            maxM = leftM;
        }

        leftN--;
        int i = kv.key;
        leftM -= g.deg[i];

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
            updateLLU<__int128>(heap, nbr, dp[nbr][h - 1]);
        }
        g.deg[i] = 0;
    }

    delete[] NotColor;
    delete[] MustColor;
}

bool tryMid(Graph &g, int h, __int128 mid, __int128 **dp, int *delArr, bool *changeFlag, int **CC, int &leftN, int *color)
{
    int delNum, nbrID;

    while (leftN)
    {
        delNum = 0;

        for (int i = 0; i < g.n; i++)
        {
            if (g.deg[i] > 0)
            {
                if (dp[i][h - 1] < mid)
                {
                    delArr[delNum++] = i;
                    for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
                    {
                        nbrID = g.adj[j];
                        if (dp[nbrID][h - 1] >= mid)
                            changeFlag[nbrID] = true;
                    }
                }
            }
        }

        if (delNum == 0)
        {
            return true;
        }

        if (delNum == leftN)
        {
            return false;
        }

        leftN -= delNum;

        g.deleteNodes(delArr, delNum);

        int cnt = 0, less = 0;
        for (int i = 0; i < g.n; i++)
        {
            if (g.deg[i] > 0)
            {
                cnt++;
                if (g.deg[i] < 113)
                    less++;
            }
        }

        for (int i = 0; i < g.n; i++)
        {
            if (changeFlag[i] == true)
            {
                changeFlag[i] = false;
                if (g.deg[i] == 0)
                {
                    leftN--; // g.deg[i] is reduced to 0 due to g.deleteNodes(delArr, delNum);
                    continue;
                }
                int colorNum_i = 0;
                for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
                {
                    nbrID = g.adj[j];
                    CC[i][color[nbrID]]++;
                    colorNum_i = max(colorNum_i, color[nbrID]);
                }

                dp[i][0] = 1;
                memset(&dp[i][1], 0, (h - 1) * sizeof(__int128));
                for (int c = 0; c <= colorNum_i; c++)
                {
                    for (int j = h - 1; j > 0; j--)
                        dp[i][j] = dp[i][j] + dp[i][j - 1] * CC[i][c];

                    CC[i][c] = 0;
                }
            }
        }
    }
    return false;
}

__int128 BinaryMaxCore(Graph &g, int h, __int128 lowerBound, __int128 upperBound, int **CC, __int128 **dp, int *color)
{

    int *degBak = new int[g.n];
    memcpy(degBak, g.deg, g.n * sizeof(int));

    int *delArr = new int[g.n];
    bool *changeFlag = new bool[g.n]();
    __int128 *dpBak = new __int128[g.n];
    int deg0 = 0;
    for (int i = 0; i < g.n; i++)
    {
        dpBak[i] = dp[i][h - 1];
        if (g.deg[i] == 0)
            deg0++;
    }

    int leftN = g.n - deg0, leftNbak = g.n;

    int times = 10;
    __int128 mid, tryValue = 100000;

    while (true)
    {
        bool result = tryMid(g, h, tryValue, dp, delArr, changeFlag, CC, leftN, color);

        if (result == false)
        {
            leftN = leftNbak;
            upperBound = mid;

            for (int i = 0; i < g.n; i++)
            {
                if (degBak[i] > 0)
                {
                    g.deg[i] = degBak[i];
                    dp[i][h - 1] = dpBak[i];
                }
            }
            upperBound = tryValue;
            break;
        }
        else
        {
            leftNbak = leftN;
            lowerBound = mid;
            for (int i = 0; i < g.n; i++)
            {
                if (degBak[i] > 0)
                {
                    degBak[i] = g.deg[i];
                    dpBak[i] = dp[i][h - 1];
                }
            }
            lowerBound = tryValue;
            tryValue *= times;
        }
        printf("lowerBound: %s\t upperBound: %s\n", _int128_to_str(lowerBound), _int128_to_str(upperBound));
    }

    while (lowerBound + 1 < upperBound)
    {
        mid = (lowerBound + upperBound) / 2;
        bool result = tryMid(g, h, mid, dp, delArr, changeFlag, CC, leftN, color);
        if (result == false)
        {
            leftN = leftNbak;
            upperBound = mid;
            for (int i = 0; i < g.n; i++)
            {
                if (degBak[i] > 0)
                {
                    g.deg[i] = degBak[i];
                    dp[i][h - 1] = dpBak[i];
                }
            }
        }
        else
        {
            leftNbak = leftN;
            lowerBound = mid;
            for (int i = 0; i < g.n; i++)
            {
                if (degBak[i] > 0)
                {
                    degBak[i] = g.deg[i];
                    dpBak[i] = dp[i][h - 1];
                }
            }
        }
        printf("mid: %s\tresult: %d\tleftN: %d\n", _int128_to_str(mid), result, leftN);
    }
    printf("leftN: %d\tlowerBound: %s\n", leftN, _int128_to_str(lowerBound));
    return lowerBound;
}

bool ubCMP(const pair<int, __int128> &a, const pair<int, __int128> &b)
{
    if (a.second == b.second)
        return a.first > b.first;
    return a.second > b.second;
}

bool ubCMPRe(const pair<int, __int128> &a, const pair<int, __int128> &b)
{
    // return a.first < b.first;
    if (a.second == b.second)
        return a.first < b.first;
    return a.second < b.second;
}

void insertSort(pair<int, __int128> *idUB, int n)
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

void compHIndex(Graph &g, int u, __int128 **dp, __int128 *uBound, int h, int *color, int **CC, int colorNum,
                pair<int, __int128> *idUB, __int128 *NotColor, __int128 *MustColor, int *lastPos)
{
    for (int i = g.cd[u]; i < g.cd[u] + g.deg[u]; i++)
    {
        int nbrID = g.adj[i];
        idUB[i - g.cd[u]].first = nbrID;
        idUB[i - g.cd[u]].second = uBound[nbrID];
    }

    sort(idUB, idUB + g.deg[u], ubCMP);

    for (int i = g.cd[u]; i < g.cd[u] + g.deg[u]; i++)
    {
        g.adj[i] = idUB[i - g.cd[u]].first;
    }

    memset(CC[u], 0, sizeof(int) * colorNum);
    memset(dp[u] + 1, 0, sizeof(__int128) * (h - 1));

    for (int L = 1; L <= g.deg[u]; L++)
    {
        int nbr = idUB[L - 1].first;
        CC[u][color[nbr]]++;

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

        if (L == g.deg[u] || !(dp[u][h - 1] < idUB[L - 1].second && dp[u][h - 1] < idUB[L].second))
        {
            uBound[u] = min(dp[u][h - 1], idUB[L - 1].second);
            lastPos[u] = L - 1;
            break;
        }
    }
}

void ColorfulStarHIndex(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 *ColofulStarCoreNum, int colorNum)
{
    bool updateFlag = true;
    // bool *flag = new bool[g.n], *lastFlag = new bool[g.n];

    __int128 *uBound = new __int128[g.n];
    for (int i = 0; i < g.n; i++)
        uBound[i] = dp[i][h - 1]; //, flag[i] = true, lastFlag[i] = true;

    // pair<int, __int128> *idUB = new pair<int, __int128>[g.maxDeg];

    pair<int, __int128> *idUB = new pair<int, __int128>[g.maxDeg];

    __int128 *NotColor = new __int128[h]();
    NotColor[0] = 1;
    __int128 *MustColor = new __int128[h]();

    int *lastPos = new int[g.n]();
    int times = 0, invokeNum = 0, sameNum = 0;

    fill(lastPos, lastPos + g.n, -1);

    while (updateFlag)
    {
        printf("updating %d\n", times++);
        updateFlag = false;

        invokeNum = 0, sameNum = 0;

        for (int i = 0; i < g.n; i++)
        {
            // int i = idBounds[ii].first;

            if (uBound[i] == 0)
            {
                continue;
            }

            if (lastPos[i] != -1)
            {
                bool updatedi = false;
                for (int j = g.cd[i]; j <= g.cd[i] + lastPos[i]; j++)
                {
                    int nbrID = g.adj[j];
                    if (uBound[nbrID] < uBound[i])
                    {
                        updatedi = true;
                        break;
                    }
                }

                if (updatedi == false)
                {
                    continue;
                }
            }

            __int128 ubBak = uBound[i];
            compHIndex(g, i, dp, uBound, h, color, CC, colorNum,
                       idUB, NotColor, MustColor, lastPos);

            invokeNum++;

            if (ubBak != uBound[i])
            {
                updateFlag = true;
            }
            else
            {
                sameNum++;
            }
        }

        printf("invoke: %d\tsame: %d\tchange: %d\n", invokeNum, sameNum, invokeNum - sameNum);
    }

    for (int i = 0; i < g.n; i++)
    {
        ColofulStarCoreNum[i] = uBound[i];
    }
}

__int128 ColorfulStarHIndexSingel(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 *ColofulStarCoreNum, int colorNum, int nodeId)
{

    static pair<int, __int128> *idUB = new pair<int, __int128>[g.maxDeg];

    static __int128 *NotColor = new __int128[h]();
    NotColor[0] = 1;
    static __int128 *MustColor = new __int128[h]();

    static int *lastPos = new int[g.n]();
    // int* lastId = new int[g.n]();

    // fill(lastId, lastId + g.n, -1);
    //__int128* boundCore = new __int128[g.n]();

    // bool prune = false;

    compHIndex(g, nodeId, dp, ColofulStarCoreNum, h, color, CC, colorNum, idUB, NotColor, MustColor, lastPos);
    //__int128 imHindex = compHIndex(g, nodeId, dp, uBound, h, nbrUB, color, CC, colorNum, idUB, NotColor, MustColor, lastPos, lastId, boundCore, &prune);

    //cout << nodeId << " -- ---- now ---- -- " << _int128_to_str(ColofulStarCoreNum[nodeId]) << endl;
    return ColofulStarCoreNum[nodeId];
}

void ColorfulStarHIndexNodes(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 *ColofulStarCoreNum, int colorNum, int numNbrBwn, int *affectedNodes)
{
    bool updateFlag = true;

    // __int128 *uBound = new __int128[g.n];
    // for (int i = 0; i < g.n; i++)
    //     uBound[i] = dp[i][h - 1];

    static pair<int, __int128> *idUB = new pair<int, __int128>[g.maxDeg];

    static __int128 *NotColor = new __int128[h]();
    NotColor[0] = 1;
    static __int128 *MustColor = new __int128[h]();

    static int *lastPos = new int[g.n]();
    for (int i = 0; i < numNbrBwn; i++)
    {
        lastPos[affectedNodes[i]] = -1;
    }

    int times = 0, invokeNum = 0, sameNum = 0;

    // fill(lastPos, lastPos + g.n, -1);

    while (updateFlag)
    {
        //printf("updating %d\n", times++);
        updateFlag = false;

        invokeNum = 0, sameNum = 0;

        for (int p = 0; p < numNbrBwn; p++)
        {
            int i = affectedNodes[p];
            // int i = idBounds[ii].first;

            if (ColofulStarCoreNum[i] == 0)
            {
                continue;
            }

            if (lastPos[i] != -1)
            {
                bool updatedi = false;
                for (int j = g.cd[i]; j <= g.cd[i] + lastPos[i]; j++)
                {
                    int nbrID = g.adj[j];
                    if (ColofulStarCoreNum[nbrID] < ColofulStarCoreNum[i])
                    {
                        updatedi = true;
                        break;
                    }
                }

                if (updatedi == false)
                {
                    continue;
                }
            }

            __int128 ubBak = ColofulStarCoreNum[i];
            compHIndex(g, i, dp, ColofulStarCoreNum, h, color, CC, colorNum,
                       idUB, NotColor, MustColor, lastPos);

            invokeNum++;

            if (ubBak != ColofulStarCoreNum[i])
            {
                updateFlag = true;
            }
            else
            {
                sameNum++;
            }
        }

        //printf("invoke: %d\tsame: %d\tchange: %d\n", invokeNum, sameNum, invokeNum - sameNum);
    }

    // for (int i = 0; i < g.n; i++)
    // {
    //     ColofulStarCoreNum[i] = uBound[i];
    // }
}

void ColStarDegreeNodes(Graph &g, __int128 **dp, int h, int colorNum, int *color, int **CC, __int128 *ColofulStarCoreNum, int numNbrBwn, int *affectedNodes, __int128 oldHindexS)
{
    for (int p = 0; p < numNbrBwn; p++)
    {
        int i = affectedNodes[p];
        dp[i][0] = 1;
        int colorNum_i = 0;
        // int* C = new int[colorNum]();

        memset(CC[i], 0, sizeof(int) * colorNum);
        memset(dp[i] + 1, 0, sizeof(__int128) * (h - 1));

        for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
        {
            int nbr = g.adj[j];
            if(ColofulStarCoreNum[nbr] < oldHindexS) continue;
            
            CC[i][color[nbr]]++;
            colorNum_i = max(colorNum_i, color[nbr]);
        }

        for (int c = 0; c <= colorNum_i; c++)
        {
            for (int j = h - 1; j > 0; j--)
                dp[i][j] = dp[i][j] + dp[i][j - 1] * CC[i][c];
        }
        ColofulStarCoreNum[i] = dp[i][h - 1];
    }
}

__int128 ColorStarDegreeWithinCore(Graph &g, __int128 **dp, int h, int colorNum, int *color, int **CC, __int128 *ColofulStarCoreNum, int sId)
{
    dp[sId][0] = 1;
    int colorNum_i = 0;
    // int* C = new int[colorNum]();

    memset(CC[sId], 0, sizeof(int) * colorNum);
    memset(dp[sId] + 1, 0, sizeof(__int128) * (h - 1));

    for (int j = g.cd[sId]; j < g.cd[sId] + g.deg[sId]; j++)
    {
        int nbr = g.adj[j];
        if(ColofulStarCoreNum[nbr] < ColofulStarCoreNum[sId]) continue;
        
        
        
        CC[sId][color[nbr]]++;
        colorNum_i = max(colorNum_i, color[nbr]);
    }

    for (int c = 0; c <= colorNum_i; c++)
    {
        for (int j = h - 1; j > 0; j--)
            dp[sId][j] = dp[sId][j] + dp[sId][j - 1] * CC[sId][c];
    }
    return dp[sId][h - 1];
    //ColofulStarCoreNum[sId] = dp[sId][h - 1];
}