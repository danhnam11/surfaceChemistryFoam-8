Foam::wordList speciesList = composition.species();
int nSpecies = speciesList.size();

//int N = mesh.nCells();
int N = numberOfMesh;
List<scalar> W(N, 1.0);
int MAXDEG = 3;
int NDEG;
double C = 0.0;
double EPS = 0.0;
int IERR;

List<scalar> R(N, 0.0);
List<scalar> muA(3 * N + 3 * MAXDEG + 3, 0.0);
List<scalar> kappaA(3 * N + 3 * MAXDEG + 3, 0.0);
List<List<scalar>> DijA(nSpecies);
forAll(DijA, i)
{
    DijA[i] = List<scalar>(3 * N + 3 * MAXDEG + 3, 0.0);
}

List<List<scalar>> muTC(nSpecies, List<scalar>(MAXDEG + 1, 0.0));
List<List<scalar>> kappaTC(nSpecies, List<scalar>(MAXDEG + 1, 0.0));
List<List<List<scalar>>> DijTC(nSpecies);
forAll(DijTC, i)
{
    DijTC[i].resize(nSpecies);
    forAll(DijTC[i], j)
    {
        DijTC[i][j].resize(MAXDEG + 1, 0.0);
    }
}

List<List<scalar>> muList(nSpecies, List<scalar>(N, 0.0));
List<List<scalar>> kappaList(nSpecies, List<scalar>(N, 0.0));
List<List<List<scalar>>> DijList(nSpecies);
forAll(DijList, i)
{
    DijList[i].resize(nSpecies);

    forAll(DijList[i], j)
    {
        DijList[i][j].resize(N, 0.0);
    }
}

List<List<scalar>> muLog = muList;
List<List<scalar>> kappaLog = kappaList;
List<List<List<scalar>>> DijLog = DijList;

List<scalar> TList(N);
List<scalar> pList(N);

List<scalar> TLog = TList;

forAll(TList, dataI)
{
    TList[dataI] = TMin + (TMax - TMin) * scalar(dataI) / scalar(N - 1);
    pList[dataI] = pressure;
}

