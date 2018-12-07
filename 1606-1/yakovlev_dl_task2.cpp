#include <mpi.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>

static const int nBaseLog = 2;

size_t TypeSize(MPI_Datatype type);
void* PointerOffset(MPI_Datatype type, void* buf, unsigned int nOffset);
int OlderBit(int nNum);
int NumberOlderBit(int nNum);
int Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

template <typename Type>
Type* CreateArray(size_t unSizeArr)
{
	Type* pArr = new Type[unSizeArr];
	for (size_t i = 0; i < unSizeArr; i++)
	{
		pArr[i] = static_cast<Type>(i);
	}
	return pArr;
}

template <>
char* CreateArray(size_t unSizeArr)
{
	char* pArr = new char[unSizeArr];
	pArr[0] = 'a';
	for (size_t i = 1; i < unSizeArr; i++)
	{
		pArr[i] = pArr[i-1] + 1;
	}
	return pArr;
}

size_t nSizeArr = 9;

int main(int argc, char**argv) 
{
	int ProcRank, ProcNum;
	int nRoot{ 0 };
	if (argc > 1)
	{
		nRoot = atoi(argv[1]);
		if (argc > 2)
			nSizeArr = atoi(argv[2]);
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (nRoot >= ProcNum)
	{
		if (ProcRank == 0)
			std::cout << "Incorect root process! Terminate!" << std::endl;
		MPI_Finalize();
		return -1;
	}
	
	char* psArr = CreateArray<char>(nSizeArr);
	int* pnArr = CreateArray<int>(nSizeArr);
	double* pdArr = CreateArray<double>(nSizeArr);
	
	char* pCharRecvBuf = new char[nSizeArr / ProcNum];
	int* pIntRecvBuf = new int[nSizeArr / ProcNum];
	double* pDoubleRecvBuf = new double[nSizeArr / ProcNum];
	
	double dStartTime{ 0 };

	if (ProcRank == nRoot)
		dStartTime = MPI_Wtime();
	MPI_Scatter(psArr, nSizeArr / ProcNum, MPI_CHAR, pCharRecvBuf, nSizeArr / ProcNum, MPI_CHAR, nRoot, MPI_COMM_WORLD);
	if (ProcRank == nRoot)
	{
		std::cout << "MPI_Scatter for datatype's int has been done in time: " << MPI_Wtime() - dStartTime << std::endl;
		dStartTime = MPI_Wtime();
	}
	MPI_Scatter(pnArr, nSizeArr / ProcNum, MPI_INT, pIntRecvBuf, nSizeArr / ProcNum, MPI_INT, nRoot, MPI_COMM_WORLD);
	if (ProcRank == nRoot)
	{
		std::cout << "MPI_Scatter for datatype's int has been done in time: " << MPI_Wtime() - dStartTime << std::endl;
		dStartTime = MPI_Wtime();
	}
	MPI_Scatter(pdArr, nSizeArr / ProcNum, MPI_DOUBLE, pDoubleRecvBuf, nSizeArr / ProcNum, MPI_DOUBLE, nRoot, MPI_COMM_WORLD);
	if (ProcRank == nRoot)
	{
		std::cout << "MPI_Scatter for datatype's int has been done in time: " << MPI_Wtime() - dStartTime << std::endl;
		dStartTime = MPI_Wtime();
	}
	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pCharRecvBuf[" << i << "] ";
		std::cout << pCharRecvBuf[i] << std::endl;
	}

	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pIntRecvBuf[" << i << "] ";
		std::cout << pIntRecvBuf[i] << std::endl;
	}

	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pDoubleRecvBuf[" << i << "] ";
		std::cout << pDoubleRecvBuf[i] << std::endl;
	}
	
	Scatter(psArr, nSizeArr / ProcNum, MPI_CHAR, pCharRecvBuf, nSizeArr / ProcNum, MPI_CHAR, nRoot, MPI_COMM_WORLD);
	if (ProcRank == nRoot)
	{
		std::cout << "MY_Scatter for datatype's char has been done in time: " << MPI_Wtime() - dStartTime << std::endl;
		dStartTime = MPI_Wtime();
	}
	
	Scatter(pnArr, nSizeArr / ProcNum, MPI_INT, pIntRecvBuf, nSizeArr / ProcNum, MPI_INT, nRoot, MPI_COMM_WORLD);
	if (ProcRank == nRoot)
	{
		std::cout << "MY_Scatter for datatype's int has been done in time: " << MPI_Wtime() - dStartTime << std::endl;
		dStartTime = MPI_Wtime();
	}
	
	Scatter(pdArr, nSizeArr / ProcNum, MPI_DOUBLE, pDoubleRecvBuf, nSizeArr / ProcNum, MPI_DOUBLE, nRoot, MPI_COMM_WORLD);
	if (ProcRank == nRoot)
	{
		std::cout << "MY_Scatter for datatype's double has been done in time: " << MPI_Wtime() - dStartTime << std::endl;
		dStartTime = MPI_Wtime();
	}
	
	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pCharRecvBuf[" << i << "] ";
		std::cout << pCharRecvBuf[i] << std::endl;
	}

	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pIntRecvBuf[" << i << "] ";
		std::cout << pIntRecvBuf[i] << std::endl;
	}
	
	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pDoubleRecvBuf[" << i << "] ";
		std::cout << pDoubleRecvBuf[i] << std::endl;
	}

	delete[] psArr, pnArr, pdArr, pCharRecvBuf, pIntRecvBuf, pDoubleRecvBuf;
	MPI_Finalize();
}


int Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
	int ProcNum, ProcRank;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	/*if (ProcRank == root)
	{
	for (size_t i = 0; i < ProcNum; i++)
	{
	if (i == root)
	{
	memcpy(recvbuf, PointerOffset(sendtype, sendbuf, sendcount*i), sendcount * TypeSize(sendtype));
	continue;
	}
	MPI_Send(PointerOffset(sendtype, sendbuf, sendcount*i), sendcount, sendtype, i, 0, comm);
	}
	if (ProcRank != root)
	MPI_Recv(recvbuf, recvcount, recvtype, root, 0, comm, &status);
	}*/

	int nRankReceiver{ 0 };
	int nRankSender{ 0 };
	int nNewProcRank = (ProcRank < root ? ProcRank + 1 : ProcRank);

	void * pvTmpBuf = sendbuf;
	
	if (ProcRank != root)
	{
		pvTmpBuf = std::malloc(sendcount * ProcNum * TypeSize(sendtype));

		int nOldBit = OlderBit(nNewProcRank);
		nRankSender = nNewProcRank & (~nOldBit);
		// те, кто раньше принимал от нулевого процесса - получают от рута
		if (nRankSender == 0)
			nRankSender = root;
		else if (nRankSender <= root)
			nRankSender--;
		MPI_Recv(pvTmpBuf, sendcount * ProcNum, recvtype, nRankSender, 0, comm, &status);
	}
	int unCicleSize = NumberOlderBit(ProcNum) - 1; // = std::log(ProcNum) / std::log(nBaseLog)
	if (ProcNum % nBaseLog > 0)
		unCicleSize++;
	assert(unCicleSize >= 0);
	int i = (ProcRank == root ? 0 : NumberOlderBit(nNewProcRank));
		
	for (; i <= unCicleSize; i++)
	{
		if (ProcRank == root)
		{
			nRankReceiver = 1 << i;
			nRankReceiver <= root ? nRankReceiver-- : nRankReceiver;
		}
		else
		{
			nRankReceiver = (1 << i) + nNewProcRank;
			if (nRankReceiver <= root)
				nRankReceiver--;
		}
		if (nRankReceiver >= ProcNum)
				break;
		assert(nRankReceiver != root);
		MPI_Send(pvTmpBuf, sendcount*ProcNum, sendtype, nRankReceiver, 0, comm);
	}
	if (ProcRank < root)
		memcpy(recvbuf, PointerOffset(sendtype, pvTmpBuf, sendcount * (ProcRank + 1)), sendcount * TypeSize(sendtype));
	else if (ProcRank == root)
		memcpy(recvbuf, sendbuf, sendcount * TypeSize(sendtype));
	else if (ProcRank > root)
		memcpy(recvbuf, PointerOffset(sendtype, pvTmpBuf, sendcount * ProcRank), sendcount * TypeSize(sendtype));

	if (ProcRank != root)
	{
		free(pvTmpBuf);
	}

	return 0;
}

size_t TypeSize(MPI_Datatype type)
{
	size_t result;
	switch (type)
	{
	default:
		std::cout << "Incorrect datatype" << std::endl;
		exit(1);
	case MPI_CHAR:
		result = sizeof(char);
		break;
	case MPI_UNSIGNED_CHAR:
		result = sizeof(unsigned char);
		break;
	case MPI_SHORT:
		result = sizeof(short);
		break;
	case MPI_UNSIGNED_SHORT:
		result = sizeof(unsigned short);
		break;
	case MPI_INT:
		result = sizeof(int);
		break;
	case MPI_UNSIGNED:
		result = sizeof(unsigned int);
		break;
	case MPI_LONG:
		result = sizeof(long);
		break;
	case MPI_UNSIGNED_LONG:
		result = sizeof(unsigned long);
		break;
	case MPI_LONG_LONG_INT:
		result = sizeof(long long);
		break;
	case MPI_FLOAT:
		result = sizeof(float);
		break;
	case MPI_DOUBLE:
		result = sizeof(double);
		break;
	case MPI_LONG_DOUBLE:
		result = sizeof(long double);
		break;
	}
	return result;
}
void* PointerOffset(MPI_Datatype type, void* buf, unsigned int nOffset)
{
	void* buffer;
	if (!buf)
	{
		std::cout << "Pointer must not be nullptr!" << std::endl;
		exit(1);
	}
	// проверить константность? 
	switch (type)
	{
	default:
		std::cout << "Incorrect datatype" << std::endl;
		exit(1);
	case MPI_CHAR:
		buffer = reinterpret_cast<char*>(buf) + nOffset;
		break;
	case MPI_UNSIGNED_CHAR:
		buffer = reinterpret_cast<unsigned char*>(buf) + nOffset;
		break;
	case MPI_SHORT:
		buffer = reinterpret_cast<short*>(buf) + nOffset;
		break;
	case MPI_UNSIGNED_SHORT:
		buffer = reinterpret_cast<unsigned short*>(buf) + nOffset;
		break;
	case MPI_INT:
		buffer = reinterpret_cast<int*>(buf) + nOffset;
		break;
	case MPI_UNSIGNED:
		buffer = reinterpret_cast<unsigned int*>(buf) + nOffset;
		break;
	case MPI_LONG:
		buffer = reinterpret_cast<long*>(buf) + nOffset;
		break;
	case MPI_UNSIGNED_LONG:
		buffer = reinterpret_cast<unsigned long*>(buf) + nOffset;
		break;
	case MPI_LONG_LONG_INT:
		buffer = reinterpret_cast<long long int*>(buf) + nOffset;
		break;
	case MPI_FLOAT:
		buffer = reinterpret_cast<float*>(buf) + nOffset;
		break;
	case MPI_DOUBLE:
		buffer = reinterpret_cast<double*>(buf) + nOffset;
		break;
	case MPI_LONG_DOUBLE:
		buffer = reinterpret_cast<long double*>(buf) + nOffset;
		break;
	}
	return buffer;
}

int OlderBit(int nNum)
{
	int nTmp = 1 << 30;
	while (nNum < nTmp)
	{
		nTmp >>= 1;
	}
	return nTmp;
}

int NumberOlderBit(int nNum)
{
	int nTmp = 1 << 30;
	int nCount = 31;
	while (nNum < nTmp)
	{
		nTmp >>= 1;
		nCount--;
	}
	return nCount;
}

