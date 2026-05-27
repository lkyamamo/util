void PairUSC::twobodyAdvanced(int i, int j, int itype, int jtype, double rsq, double &fpair, int eflag, double &eng, bool skipAtomPair) {
    double Ni0 = coord.coord[i][0];
    double Ni1 = coord.coord[i][1];

    double Nj0 = coord.coord[j][0];
    double Nj1 = coord.coord[j][1];

    if(Ni1 == 0 && Ni0 == 0 && Nj1 == 0 && Nj0 == 0) {
        // Assume first type if there are no other atoms
        Nj0 = 1.0;
        Ni0 = 1.0;
    }
    const double oneOverCoordinateNumbers = 1.0 / (Ni0 + Nj0 + Ni1 + Nj1);
    const double g = (1.0 - (Ni1 + Nj1) * oneOverCoordinateNumbers);
    dg[coord.typesInPairs[0][1]] = -(Ni1 + Nj1)*oneOverCoordinateNumbers*oneOverCoordinateNumbers;
    dg[coord.typesInPairs[1][1]] = (Ni0 + Nj0)*oneOverCoordinateNumbers*oneOverCoordinateNumbers;

    // Compute advanced interpolated force and energy
    double fpair1 = 0.0;
    double evdwl1 = 0.0;
    double fpair2 = 0.0;
    double evdwl2 = 0.0;
    int ijparam = elem2param[itype][jtype][jtype][0];
    twobody(&params[ijparam],rsq,fpair1,eflag,evdwl1);
    ijparam = elem2param[itype][jtype][jtype][1];
    twobody(&params[ijparam],rsq,fpair2,eflag,evdwl2);

    // This is the interpolation weight between the two types
    if(!skipAtomPair) {
        fpair = g*fpair1 + (1.0-g)*fpair2;
        eng = + g*evdwl1 + (1.0-g)*evdwl2;
    } else {
        fpair = 0;
        eng = 0;
    }

    if(std::min(1.0-g, g) > 1e-1) {
        twobodyGradF(i, j, skipAtomPair, itype, jtype, evdwl1, evdwl2);
    }
}