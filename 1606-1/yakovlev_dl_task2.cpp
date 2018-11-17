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


static size_t nSizeArr = 9;

int main(int argc, char**argv)
{
	char str[]{ "abcdefghi" };
	int nums[]{ 1,2,3,4,5,6,7,8,9 };
	int ProcRank, ProcNum;
	int nRoot;
	if (argc > 1)
		nRoot = atoi(argv[1]);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	char* pCharRecvBuf = new char[nSizeArr / ProcNum];
	int* pIntRecvBuf = new int[nSizeArr / ProcNum];
	MPI_Scatter(str, nSizeArr / ProcNum, MPI_CHAR, pCharRecvBuf, nSizeArr / ProcNum, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(nums, nSizeArr / ProcNum, MPI_INT, pIntRecvBuf, nSizeArr / ProcNum, MPI_INT, nRoot, MPI_COMM_WORLD);
	
	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pIntRecvBuf[" << i << "] ";
		std::cout << pIntRecvBuf[i] << std::endl;
	}
	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pCharRecvBuf[" << i << "] ";
		std::cout << pCharRecvBuf[i] << std::endl;
	}
	
	Scatter(str, nSizeArr / ProcNum, MPI_CHAR, pCharRecvBuf, nSizeArr / ProcNum, MPI_CHAR, 0, MPI_COMM_WORLD);
	Scatter(nums, nSizeArr / ProcNum, MPI_INT, pIntRecvBuf, nSizeArr / ProcNum, MPI_INT, nRoot, MPI_COMM_WORLD);
	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pIntRecvBuf[" << i << "] ";
		std::cout << pIntRecvBuf[i] << std::endl;
	}
	for (size_t i = 0; i < nSizeArr / ProcNum; i++)
	{
		std::cout << "ProcRank " << ProcRank << " pCharRecvBuf[" << i << "] ";
		std::cout << pCharRecvBuf[i] << std::endl;
	}
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
	int unCicleSize = NumberOlderBit(ProcNum) - 1; // = std::log(ProcNum) / std::log(nBaseLog);
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
