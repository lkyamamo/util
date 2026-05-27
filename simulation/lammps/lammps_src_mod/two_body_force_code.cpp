void PairUSC::twobodyGradF(int i, int j, bool skipAtomPair, int itype, int jtype, double energy1, double energy2) {
    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    // Neighbor list for atom i
    const int *neighborlist = firstneigh[i];
    const int numNeighbors = numneigh[i];

    // These types are 1 smaller since they are used as array indices earlier through the map
    double xi[3];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];

    for (int kk = 0; kk < numNeighbors; kk++) {
        const int k = neighborlist[kk] & NEIGHMASK;

        int ktype = map[type[k]];
        if(ktype == itype) continue;
        int atomPairIndex = coord.coordTypeIndexMap[itype][ktype];
        if(atomPairIndex < 0) {
            continue;
        }

        const double dx = xi[0] - x[k][0];
        const double dy = xi[1] - x[k][1];
        const double dz = xi[2] - x[k][2];
        const double rsq = dx*dx + dy*dy + dz*dz;

        // The derivative of the interpolation function theta(r)
        if(withinRange(rsq, coord.R[atomPairIndex], coord.D[atomPairIndex])) {
            int potentialTableIndex = (rsq-tabinnersq)*oneOverDeltaR2;
            float fractionalRest = (rsq-tabinnersq)*oneOverDeltaR2 - potentialTableIndex; // double - int will only keep the 0.xxxx part
            float gradThetaOverR1 = gradThetaOverRTable[atomPairIndex][potentialTableIndex];
            float gradThetaOverR2 = gradThetaOverRTable[atomPairIndex][potentialTableIndex+1];
            float gradThetaOverR = (1.0 - fractionalRest)*gradThetaOverR1 + fractionalRest*gradThetaOverR2;

            const float dgTimesGradTheta = dg[ktype]*gradThetaOverR;
            const float F = dgTimesGradTheta*(energy1 - energy2);

            f[i][0] += F*dx;
            f[i][1] += F*dy;
            f[i][2] += F*dz;
            f[k][0] -= F*dx;
            f[k][1] -= F*dy;
            f[k][2] -= F*dz;
        }
    }

    if(j < atom->nlocal && !skipAtomPair) {
        // Neighbor list for atom j
        const int *neighborlist = firstneigh[j];
        const int numNeighbors = numneigh[j];
        double xj[3];
        xj[0] = x[j][0];
        xj[1] = x[j][1];
        xj[2] = x[j][2];

        for (int kk = 0; kk < numNeighbors; kk++) {
            const int k = neighborlist[kk] & NEIGHMASK;
            int ktype = map[type[k]];
            if(ktype == jtype) continue;
            int atomPairIndex = coord.coordTypeIndexMap[jtype][ktype];
            if(atomPairIndex < 0) {
                continue;
            }

            const double dx = xj[0] - x[k][0];
            const double dy = xj[1] - x[k][1];
            const double dz = xj[2] - x[k][2];
            const double rsq = dx*dx + dy*dy + dz*dz;

            // The derivative of the interpolation function theta(r)
            if(withinRange(rsq, coord.R[atomPairIndex], coord.D[atomPairIndex])) {
                int potentialTableIndex = (rsq-tabinnersq)*oneOverDeltaR2;
                float fractionalRest = (rsq-tabinnersq)*oneOverDeltaR2 - potentialTableIndex; // double - int will only keep the 0.xxxx part
                float gradThetaOverR1 = gradThetaOverRTable[atomPairIndex][potentialTableIndex];
                float gradThetaOverR2 = gradThetaOverRTable[atomPairIndex][potentialTableIndex+1];
                float gradThetaOverR = (1.0 - fractionalRest)*gradThetaOverR1 + fractionalRest*gradThetaOverR2;

                const float dgTimesGradTheta = dg[ktype]*gradThetaOverR;
                const float F = dgTimesGradTheta*(energy1 - energy2);

                f[j][0] += F*dx;
                f[j][1] += F*dy;
                f[j][2] += F*dz;
                f[k][0] -= F*dx;
                f[k][1] -= F*dy;
                f[k][2] -= F*dz;
            }
        }
    }
}