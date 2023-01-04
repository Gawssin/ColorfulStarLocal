void initColStarDegree(Graph &g, __int128 **dp, int h, int colorNum, int *color, int **CC, int basicFlag = 1)
{
    __int128 *NotColor0 = new __int128[h]();
    __int128 *MustColor0 = new __int128[h]();

#pragma omp parallel for schedule(dynamic, 16)
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

bool ubCMP(const pair<int, __int128> &a, const pair<int, __int128> &b)
{
    if (a.second == b.second)
        return a.first > b.first;
    return a.second > b.second;
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

void ColorfulStarHIndex(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 *ColofulStarCoreNum, int colorNum, int optimization)
{
    bool updateFlag = true;

    __int128 *uBound = new __int128[g.n];
    for (int i = 0; i < g.n; i++)
        uBound[i] = dp[i][h - 1];

    int numThreads = omp_get_max_threads();
    cout << "number of Threads: " << numThreads << endl;

    pair<int, __int128> **idUBThreads = new pair<int, __int128> *[numThreads];

    __int128 **NotColorThreads = new __int128 *[numThreads];
    __int128 **MustColorThreads = new __int128 *[numThreads];

    for (int i = 0; i < numThreads; i++)
    {
        idUBThreads[i] = new pair<int, __int128>[g.maxDeg];

        NotColorThreads[i] = new __int128[h]();
        NotColorThreads[i][0] = 1;
        MustColorThreads[i] = new __int128[h]();
    }

    int *lastPos = new int[g.n]();
    int times = 0, invokeNum = 0, sameNum = 0;

    fill(lastPos, lastPos + g.n, -1);

    long long totalInvoke = 0;
    while (updateFlag)
    {
        printf("updating %d\n", times++);
        updateFlag = false;

        invokeNum = 0, sameNum = 0;

#pragma omp parallel for schedule(dynamic, 16) reduction(+ \
                                                         : invokeNum, sameNum)
        for (int i = 0; i < g.n; i++)
        {
            if (uBound[i] == 0)
            {
                continue;
            }

            if (optimization >= 3)
            {

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
            }

            int threadId = omp_get_thread_num();

            __int128 ubBak = uBound[i];
            // compute H-index for node i
            compHIndex(g, i, dp, uBound, h, color, CC, colorNum,
                       idUBThreads[threadId], NotColorThreads[threadId], MustColorThreads[threadId], lastPos);

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

        totalInvoke += invokeNum;

        printf("invoke: %d\tsame: %d\tchange: %d\n", invokeNum, sameNum, invokeNum - sameNum);
    }

    printf("Total invoke:\t %lld\n", totalInvoke);
    printf("Total times (iterations):\t %d\n", times);

    for (int i = 0; i < g.n; i++)
    {
        ColofulStarCoreNum[i] = uBound[i];
    }
}

void compHIndexSync(Graph &g, int u, __int128 **dp, __int128 *uBound, int h, int *color, int **CC, int colorNum,
                    pair<int, __int128> *idUB, __int128 *NotColor, __int128 *MustColor, __int128 *uBoundLast, int *lastPos)
{
    for (int i = g.cd[u]; i < g.cd[u] + g.deg[u]; i++)
    {
        int nbrID = g.adj[i];
        idUB[i - g.cd[u]].first = nbrID;
        idUB[i - g.cd[u]].second = uBoundLast[nbrID];
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

void ColorfulStarHIndexSync(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 *ColofulStarCoreNum, int colorNum)
{
    bool updateFlag = true;

    __int128 *uBoundLast = new __int128[g.n];
    __int128 *uBound = new __int128[g.n];
    for (int i = 0; i < g.n; i++)
        uBoundLast[i] = dp[i][h - 1];

    int numThreads = omp_get_max_threads();
    cout << "number of Threads: " << numThreads << endl;

    pair<int, __int128> **idUBThreads = new pair<int, __int128> *[numThreads];

    __int128 **NotColorThreads = new __int128 *[numThreads];
    __int128 **MustColorThreads = new __int128 *[numThreads];

    for (int i = 0; i < numThreads; i++)
    {
        idUBThreads[i] = new pair<int, __int128>[g.maxDeg];

        NotColorThreads[i] = new __int128[h]();
        NotColorThreads[i][0] = 1;
        MustColorThreads[i] = new __int128[h]();
    }

    int times = 0, invokeNum = 0, sameNum = 0;

    int *lastPos = new int[g.n]();
    fill(lastPos, lastPos + g.n, -1);

    long long totalInvoke = 0;

    while (updateFlag)
    {
        printf("updating %d\n", times++);
        updateFlag = false;

        invokeNum = 0, sameNum = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ \
                                                         : invokeNum, sameNum)
        for (int i = 0; i < g.n; i++)
        {

            if (uBoundLast[i] == 0)
            {
                uBound[i] = 0;
                continue;
            }

            int threadId = omp_get_thread_num();

            compHIndexSync(g, i, dp, uBound, h, color, CC, colorNum, idUBThreads[threadId], NotColorThreads[threadId], MustColorThreads[threadId],
                           uBoundLast, lastPos);

            invokeNum++;

            if (uBoundLast[i] != uBound[i])
            {
                updateFlag = true;
            }
            else
            {
                sameNum++;
            }
        }

        swap(uBound, uBoundLast);
        totalInvoke += invokeNum;

        printf("invoke: %d\tsame: %d\tchange: %d\n", invokeNum, sameNum, invokeNum - sameNum);
    }

    printf("Total invoke:\t %lld\n", totalInvoke);
    printf("Total times:\t %d\n", times);

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
    compHIndex(g, nodeId, dp, ColofulStarCoreNum, h, color, CC, colorNum, idUB, NotColor, MustColor, lastPos);
    return ColofulStarCoreNum[nodeId];
}

void ColorfulStarHIndexNodes(Graph &g, __int128 **dp, int h, int *color, int **CC, __int128 *ColofulStarCoreNum, int colorNum, int numNbrBwn, int *affectedNodes)
{
    bool updateFlag = true;

    int numThreads = omp_get_max_threads();

    static pair<int, __int128> **idUBThreads = new pair<int, __int128> *[numThreads]();

    static __int128 **NotColorThreads = new __int128 *[numThreads]();
    static __int128 **MustColorThreads = new __int128 *[numThreads]();

    for (int i = 0; i < numThreads; i++)
    {
        if (idUBThreads[i] == NULL)
        {
            idUBThreads[i] = new pair<int, __int128>[g.maxDeg];
            NotColorThreads[i] = new __int128[h]();
            MustColorThreads[i] = new __int128[h]();
        }
        NotColorThreads[i][0] = 1;
    }

    static int *lastPos = new int[g.n]();
    for (int i = 0; i < numNbrBwn; i++)
    {
        lastPos[affectedNodes[i]] = -1;
    }

    int times = 0;

    while (updateFlag)
    {
        updateFlag = false;
        
        #pragma omp parallel for schedule(dynamic, 16) 
        for (int p = 0; p < numNbrBwn; p++)
        {
            int i = affectedNodes[p];
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
            
            int threadId = omp_get_thread_num();
            
            __int128 ubBak = ColofulStarCoreNum[i];
            compHIndex(g, i, dp, ColofulStarCoreNum, h, color, CC, colorNum,
                       idUBThreads[threadId], NotColorThreads[threadId], MustColorThreads[threadId], lastPos);

            if (ubBak != ColofulStarCoreNum[i])
            {
                updateFlag = true;
            }
        }
    }
}

void ColStarDegreeNodes(Graph &g, __int128 **dp, int h, int colorNum, int *color, int **CC, __int128 *ColofulStarCoreNum, int numNbrBwn,
                        int *affectedNodes, __int128 oldHindexS)
{
    #pragma omp parallel for schedule(dynamic, 16) 
    for (int p = 0; p < numNbrBwn; p++)
    {
        int i = affectedNodes[p];
        dp[i][0] = 1;
        int colorNum_i = 0;

        memset(CC[i], 0, sizeof(int) * colorNum);
        memset(dp[i] + 1, 0, sizeof(__int128) * (h - 1));

        for (int j = g.cd[i]; j < g.cd[i] + g.deg[i]; j++)
        {
            int nbr = g.adj[j];
            if (ColofulStarCoreNum[nbr] < oldHindexS)
                continue;

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

    memset(CC[sId], 0, sizeof(int) * colorNum);
    memset(dp[sId] + 1, 0, sizeof(__int128) * (h - 1));

    for (int j = g.cd[sId]; j < g.cd[sId] + g.deg[sId]; j++)
    {
        int nbr = g.adj[j];
        if (ColofulStarCoreNum[nbr] < ColofulStarCoreNum[sId])
            continue;

        CC[sId][color[nbr]]++;
        colorNum_i = max(colorNum_i, color[nbr]);
    }

    for (int c = 0; c <= colorNum_i; c++)
    {
        for (int j = h - 1; j > 0; j--)
            dp[sId][j] = dp[sId][j] + dp[sId][j - 1] * CC[sId][c];
    }
    return dp[sId][h - 1];
}