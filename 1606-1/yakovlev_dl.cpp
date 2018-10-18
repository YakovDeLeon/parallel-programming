#include <mpi.h>
#include <iostream>
#include <memory>
#include <ctime>

void FillArray(double * pd2dArray, size_t unRows, size_t unCols);
void ShowArray(double * pd2dArray, size_t unRows, size_t unCols);

int main(int argc, char** argv)
{
	int nRows = 10;
	int nCols = 30;
	
	if (argc == 3)
	{
		nRows = atoi(argv[1]);
		nCols = atoi(argv[2]);
	}
	srand(time(nullptr));
	int ProcRank, ProcNum;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	double * pdMatrix, *pdRecvMatrix, *pdChunksOfSumArray;
	double * pdSumArray = new double[nCols] {0};
	double dStartTime;
	if (ProcRank == 0)
	{
		pdMatrix = new double[nCols * nRows];
		FillArray(pdMatrix, nRows, nCols);
		ShowArray(pdMatrix, nRows, nCols);
		dStartTime = MPI_Wtime();
		std::cout << "Computer has been beginning a solving!" << std::endl;
	}
	else
	{
		pdMatrix = nullptr;
	}


	int nRemColsPerProc = nCols % ProcNum;
	
	int * pnSendCounts = new int[ProcNum];
	int * pnRecvCounts = new int[ProcNum];
	int * pnDispls = new int[ProcNum];
	int * pnRecvDispls = new int[ProcNum];
	int nSendSum = 0;
	int nRecvSum = 0;

	for (size_t i = 0; i < ProcNum; i++)
	{
		pnSendCounts[i] = nCols / ProcNum * nRows ;
		pnRecvCounts[i] = pnSendCounts[i] / nRows;
		if (nRemColsPerProc > 0)
		{
			pnRecvCounts[i]++;
			pnSendCounts[i] += nRows;
			nRemColsPerProc--;
		}
		pnDispls[i] = nSendSum;
		nSendSum += pnSendCounts[i];
		pnRecvDispls[i] = nRecvSum;
		nRecvSum += pnRecvCounts[i];
	}
	pdChunksOfSumArray = new double[pnRecvCounts[ProcRank]]{0};
	pdRecvMatrix = new double[pnSendCounts[ProcRank]];
	
	MPI_Scatterv(pdMatrix, pnSendCounts, pnDispls, MPI_DOUBLE, pdRecvMatrix, pnSendCounts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int nCount = 0;
	for (int i = 0; i < pnSendCounts[ProcRank]; i++)
	{
		if (i !=0 && i % nRows == 0)
			nCount++;
		pdChunksOfSumArray[nCount] += pdRecvMatrix[i];
	}
	MPI_Gatherv(pdChunksOfSumArray, nCount + 1, MPI_DOUBLE, pdSumArray, pnRecvCounts, pnRecvDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		double dTimeDiff = MPI_Wtime() - dStartTime;
		std::cout << "Time is " << dTimeDiff << " sec" << std::endl;
	}
	MPI_Finalize();
	
	if (ProcRank == 0)
	{
		for (size_t i = 0; i < nCols; i++)
		{
			std::cout << pdSumArray[i] << " ";
		}
		std::cout << std::endl;
	}
	delete[] pdMatrix, pdRecvMatrix, pdChunksOfSumArray, pdSumArray;
	delete[] pnSendCounts, pnRecvCounts, pnDispls, pnRecvDispls;
	return 0;
}

void FillArray(double * pd2dArray, size_t unRows, size_t unCols)
{
	for (size_t i = 0; i < unRows * unCols; ++i)
	{
			pd2dArray[i] = rand() % 100;
	}
}

void ShowArray(double * pd2dArray, size_t unRows, size_t unCols)
{
	for (size_t i = 0; i < unRows; ++i)
	{
		for (size_t j = 0; j < unCols * unRows; j += unRows)
			std::cout << pd2dArray[i + j] << " ";
		std::cout << std::endl;
	}
}
