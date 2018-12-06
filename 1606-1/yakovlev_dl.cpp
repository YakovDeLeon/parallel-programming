#include <mpi.h>
#include <iostream>
#include <memory>
#include <ctime>
#include <limits>


#pragma warning (disable: 4703)

void FillArray(double * pd2dArray, size_t unRows, size_t unCols, double dBegin, double dEnd);
void ShowArray(double * pd2dArray, size_t unRows, size_t unCols);


int main(int argc, char** argv)
{
	int nRows, nCols;
	int ProcRank, ProcNum;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0)
	{
		if (argc == 5)
		{
			nRows = atoi(argv[1]);
			nCols = atoi(argv[2]);
		}
		else
		{
			nRows = 3;
			nCols = 3;
		}
	}

	MPI_Bcast(&nRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	MPI_Datatype MatrixCol;
	MPI_Type_vector(nRows, 1, 1, MPI_DOUBLE, &MatrixCol);
	MPI_Type_commit(&MatrixCol);

	double * pdMatrix, *pdRecvMatrix, *pdChunksOfSumArray;
	double * pdSumArray;
	
	if (nCols >= ProcNum || (ProcRank == 0 && nCols < ProcNum))
		pdSumArray = new double[nCols] {0};


	double dStartTime;
	if (ProcRank == 0)
	{
		pdMatrix = new double[nCols * nRows];
		if (argc == 5)
			FillArray(pdMatrix, nRows, nCols, atof(argv[3]), atof(argv[4]));
		else
			FillArray(pdMatrix, nRows, nCols, 0.f, 10.f);
		ShowArray(pdMatrix, nRows, nCols);
		dStartTime = MPI_Wtime();
		std::cout << "Computer has been beginning a solving!" << std::endl;
	}
	else
	{
		pdMatrix = nullptr;
	}

	if (nCols < ProcNum && ProcRank == 0)
	{
		int nCount = 0;
		for (int i = 0; i < nCols * nRows; i++)
		{
			if (i != 0 && i % nRows == 0)
				nCount++;
			pdSumArray[nCount] += pdMatrix[i];
		}

	}
	else
	{
		int * pnColCounts = new int[ProcNum];
		int * pnColDispls = new int[ProcNum];
		if (ProcRank == 0)
		{
			int nSendSum = 0;
			int nRecvSum = 0;
			int nRemColsPerProc = nCols % ProcNum;
			for (size_t i = 0; i < ProcNum; i++)
			{
				pnColCounts[i] = nCols / ProcNum;
				if (nRemColsPerProc > 0)
				{
					pnColCounts[i]++;
					nRemColsPerProc--;
				}
				pnColDispls[i] = nRecvSum;
				nRecvSum += pnColCounts[i];
			}
		}
		MPI_Bcast(pnColCounts, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(pnColDispls, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);

		pdChunksOfSumArray = new double[pnColCounts[ProcRank]]{ 0 };
		pdRecvMatrix = new double[pnColCounts[ProcRank] * nRows];

		MPI_Scatterv(pdMatrix, pnColCounts, pnColDispls, MatrixCol, pdRecvMatrix, pnColCounts[ProcRank], MatrixCol, 0, MPI_COMM_WORLD);

		int nCount = 0;
		for (int i = 0; i < pnColCounts[ProcRank] * nRows; i++)
		{
			if (i != 0 && i % nRows == 0)
				nCount++;
			pdChunksOfSumArray[nCount] += pdRecvMatrix[i];
		}

		MPI_Gatherv(pdChunksOfSumArray, nCount + 1, MPI_DOUBLE, pdSumArray, pnColCounts, pnColDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		delete[] pnColCounts, pnColDispls;
	}
	if (ProcRank == 0)
	{
		double dTimeDiff = MPI_Wtime() - dStartTime;
		std::cout << "Time is " << dTimeDiff << " sec" << std::endl;
	}

	MPI_Finalize();
	if (ProcRank == 0)
	{
		ShowArray(pdSumArray, 1, nCols);
		std::cout << std::endl;
	}
	delete[] pdMatrix, pdRecvMatrix, pdChunksOfSumArray, pdSumArray;
	
	return 0;
}

void FillArray(double * pd2dArray, size_t unRows, size_t unCols, double dBegin, double dEnd)
{
	srand(time(nullptr));
	for (size_t i = 0; i < unRows * unCols; ++i)
	{
		pd2dArray[i] = rand() % (int)(dEnd - dBegin) + dBegin + double(rand()) / 1000000.f;
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
