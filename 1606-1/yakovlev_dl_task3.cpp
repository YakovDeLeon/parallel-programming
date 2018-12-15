#include <mpi.h>
#include <iostream>
#include <memory>
#include <ctime>
#include <limits>
#include <cassert>

#pragma warning (disable: 4703)

void FillArray(double * pd2dArray, size_t unSizeArr, double dBegin, double dEnd);
void ShowArray(double * pdArray, size_t unSizeArr, int ProcRank);
void* PointerOffset(MPI_Datatype type, void* buf, unsigned int nOffset);
int OlderBit(int nNum);
int NumberOlderBit(int nNum);
template <typename Type>
int Shell_Reduce(Type *sendbuf, int *sendcount, Type *recvbuf, int sizearr, MPI_Datatype type, int root, MPI_Comm comm);


int main(int argc, char** argv)
{
	int nSizeArr;
	int ProcRank, ProcNum;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0)
	{
		if (argc == 4)
		{
			nSizeArr = atoi(argv[1]);
		}
		else
		{
			nSizeArr = 1000;
		}
	}

	MPI_Bcast(&nSizeArr, 1, MPI_INT, 0, MPI_COMM_WORLD);

	double * pdArr, *pdRecvArr;
	double * pdRes = new double[nSizeArr] {0};

	/*if (nSizeArr >= ProcNum || (ProcRank == 0 && nCols < ProcNum))
		pdSumArray = new double[nCols] {0};
	*/

	double dStartTime;
	if (ProcRank == 0)
	{
		pdArr = new double[nSizeArr];
		if (argc == 5)
			FillArray(pdArr, nSizeArr, atof(argv[3]), atof(argv[4]));
		else
			FillArray(pdArr, nSizeArr, 0.f, 10.f);
		ShowArray(pdArr, nSizeArr, 0);
		dStartTime = MPI_Wtime();
		std::cout << "Computer has been beginning a solving!" << std::endl;
	}
	else
	{
		pdArr = nullptr;
	}

	/*if (nCols < ProcNum && ProcRank == 0)
	{
		int nCount = 0;
		for (int i = 0; i < nCols * nRows; i++)
		{
			if (i != 0 && i % nRows == 0)
				nCount++;
			pdSumArray[nCount] += pdMatrix[i];
		}

	}*/

	
	int * pnArrCounts = new int[ProcNum];
	int * pnArrDispls = new int[ProcNum];
	if (ProcRank == 0)
	{
		int nSendSum = 0;
		int nRecvSum = 0;
		int nRemSizePerProc = nSizeArr % ProcNum;
		for (size_t i = 0; i < ProcNum; i++)
		{
			pnArrCounts[i] = nSizeArr / ProcNum ;
			if (nRemSizePerProc > 0)
			{
				pnArrCounts[i]++;
				nRemSizePerProc--;
			}
			pnArrDispls[i] = nRecvSum;
			nRecvSum += pnArrCounts[i];
		}
	}
	MPI_Bcast(pnArrCounts, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(pnArrDispls, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);

	//pdChunksOfSumArray = new double[pnColCounts[ProcRank]]{ 0 };
	pdRecvArr = new double[pnArrCounts[ProcRank]];

	MPI_Scatterv(pdArr, pnArrCounts, pnArrDispls, MPI_DOUBLE, pdRecvArr, pnArrCounts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*int nCount = 0;
	for (int i = 0; i < pnColCounts[ProcRank] * nRows; i++)
	{
		if (i != 0 && i % nRows == 0)
			nCount++;
		pdChunksOfSumArray[nCount] += pdRecvMatrix[i];
	}

	MPI_Gatherv(pdChunksOfSumArray, nCount + 1, MPI_DOUBLE, pdSumArray, pnColCounts, pnColDispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	*/
	Shell_Reduce(pdRecvArr, pnArrCounts, pdRes, nSizeArr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete[] pnArrCounts, pnArrDispls, pdRecvArr;
	
	if (ProcRank == 0)
	{
		double dTimeDiff = MPI_Wtime() - dStartTime;
		std::cout << "Time is " << dTimeDiff << " sec" << std::endl;
	}

	MPI_Finalize();
	if (ProcRank == 0)
	{
		ShowArray(pdRes, nSizeArr, 0);
		std::cout << std::endl;
	}
	delete[] pdArr, pdRes;
	return 0;
}

void FillArray(double * pd2dArray, size_t unSizeArr, double dBegin, double dEnd)
{
	srand(time(nullptr));
	for (size_t i = 0; i < unSizeArr; ++i)
	{
		pd2dArray[i] = rand() % (int)(dEnd - dBegin) + dBegin + double(rand()) / 1000000.f;
	}
}

void ShowArray(double * pdArray, size_t unSizeArr, int ProcRank)
{
	for (size_t i = 0; i < unSizeArr; ++i)
	{
		std::cout << "ProcRank " << ProcRank << " has a folowing arr[" <<  i << "] = " << pdArray[i] << std::endl;
	}
}


template<typename Type>
Type* RealocArray(Type* pArr, size_t nOldSize, size_t nNewSize)
{
	assert(nOldSize <= nNewSize);
	Type* pNewArr = new Type[nNewSize]{ 0 };
	for (size_t i = 0; i < nOldSize; i++)
	{
		pNewArr[i] = pArr[i];
	}
	delete[] pArr;
	pArr = nullptr;
	return pNewArr;
}

template <typename BaseType>
void ShellsSort(BaseType *A, unsigned N)
{
	unsigned i, j, k;
	BaseType t;
	for (k = N / 2; k > 0; k /= 2)
		for (i = k; i < N; i++)
		{
			t = A[i];
			for (j = i; j >= k; j -= k)
			{
				if (t < A[j - k])
					A[j] = A[j - k];
				else
					break;
			}
			A[j] = t;
		}
}

/*template<typename T>
bool compare(const T a, const T b)
{
return a < b;
}
template<>
bool compare<const char*>(const char * a, const char * b)
{
if (strlen(a) == strlen(b))
return strcmp(a, b) < 0;
else
return  strlen(a) < strlen(b);
}

template<>
bool compare<char*>(char * a, char * b)
{
if (strlen(a) == strlen(b))
return strcmp(a, b) < 0;
else
return  strlen(a) < strlen(b);
}

template<typename T>
void copy(T * begin, T * end, T * dst)
{
for (size_t i = 0; begin + i != end; ++i)
*(dst + i) = *(begin + i);
}
*/

template<typename T>
void Merge(T* mas1, T* mas2, T* tmp, size_t size1, size_t size2)
{
	int a = 0;
	int b = 0;
	int i = 0;
	while ((a != size1) && (b != size2))
	{
		if (mas1[a] <= mas2[b])
		{
			tmp[i] = mas1[a];
			a++;
		}
		else
		{
			tmp[i] = mas2[b];
			b++;
		}
		i++;
	}
	if (a == size1)
	{
		int j = b;
		for (; j<size2; j++, i++)
			tmp[i] = mas2[j];
	}
	else
	{
		int j = a;
		for (; j<size1; j++, i++)
			tmp[i] = mas1[j];
	}
}


template<typename T>
void MergeArrays(T* mas1, T* mas2, T* tmp, size_t size1, size_t size2)
{
	int nMedIndex = BinSearch(mas2, 0, size2, (mas1[size1 / 2]));
	//std::cout << "Medium index = " << nMedIndex << std::endl;
	int nTmp1Size = nMedIndex + size1 / 2;
	int nTmp2Size = size1 + size2 - nTmp1Size;
	T* tmp1 = new T[nTmp1Size];
	T* tmp2 = new T[nTmp2Size];
	Merge(mas1, mas2, tmp1, size1 / 2, nMedIndex);
	Merge(mas1 + size1 / 2, mas2 + nMedIndex , tmp2, size1 - size1 / 2, size2 - (nMedIndex));
	memcpy(tmp, tmp1, nTmp1Size * sizeof(T));
	memcpy(tmp + nTmp1Size, tmp2, nTmp2Size * sizeof(T));
	delete[] tmp1, tmp2;
}

template<typename T>
int BinSearch(T *mas, int l, int r, T x)
{
	if (l == r)
		return l;
	if (l + 1 == r)
		if (x < mas[l])
			return l;
		else
			return r;
	int m = (l + r) / 2;
	if (x < mas[m])
		r = m;
	else
		if (x > mas[m])
			l = m;
		else
			return m;
	return BinSearch(mas, l, r, x);
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

template <typename Type>
int Shell_Reduce(Type *sendbuf, int *sendcount, Type *recvbuf, int sizearr, MPI_Datatype type, int root, MPI_Comm comm)
{
	MPI_Status status;
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int nRankSender{ 0 };
	Type * pvResBuf = new Type[sizearr];
	//recvbuf = malloc(sendcount[ProcRank] * TypeSize(type));
	int nFinish = std::log(ProcNum) / std::log(2);
	if (ProcNum % 2 > 0)
		nFinish++;
	ShellsSort((double*)sendbuf, sendcount[ProcRank]);
	ShowArray(sendbuf, sendcount[ProcRank], ProcRank);
	int newProcRank = (ProcRank < root ? ProcRank + 1 : ProcRank);
	int nStart = (ProcRank == root ? 0 : NumberOlderBit(newProcRank));
	for (int i = nFinish; i >= nStart; i--)
	{
		if (ProcRank == root)
		{
			nRankSender = std::pow(2, i);
			nRankSender <= root ? nRankSender-- : nRankSender;
		}
		else
		{
			nRankSender = std::pow(2, i) + newProcRank;
			if (nRankSender <= root)
				nRankSender--;
		}
		if (nRankSender < ProcNum)
		{
			sendbuf = RealocArray(sendbuf, sendcount[ProcRank], sendcount[ProcRank] + sendcount[nRankSender]);
			sendcount[ProcRank] += sendcount[nRankSender];
			std::cout << ProcRank << " recv from " << nRankSender << " array with size " << sendcount[ProcRank] <<std::endl;
			MPI_Recv(recvbuf, sendcount[ProcRank], type, nRankSender, 0, comm, &status);
			MergeArrays(sendbuf, recvbuf, pvResBuf, sendcount[ProcRank] - sendcount[nRankSender], sendcount[nRankSender]);
			
			//recvbuf = realloc(sendbuf, sendcount[ProcRank] * TypeSize(type));
			memcpy(sendbuf, pvResBuf, (sendcount[ProcRank]) * TypeSize(type));
			std::cout << "	ProcRank " << ProcRank << " has a following array:" << std::endl;
			ShowArray(sendbuf, sendcount[ProcRank], ProcRank);
			std::cout << std::endl;
			memset(recvbuf, 0, sendcount[ProcRank] * TypeSize(type));
			//Operation(type, op, sendbuf, recvbuf, sendbuf, count);
		}
	}
	if (ProcRank != root)
	{
		int oldnum = OlderBit(newProcRank);
		int nRankReciever = newProcRank & (~oldnum);
		if (nRankReciever == 0)
			nRankReciever = root;
		else if (nRankReciever <= root)
			nRankReciever--;
		std::cout << ProcRank << " send to " << nRankReciever << " array with size " << sendcount[ProcRank] << std::endl;
		MPI_Send(sendbuf, sendcount[ProcRank], type, nRankReciever, 0, comm);
	}
	if (ProcRank == root)
	{
		memcpy(recvbuf, sendbuf, sizearr*TypeSize(type));
	}

	//ShowArray((double*)recvbuf, sendcount[ProcRank]);
	//free(pvResBuf);
	delete[] pvResBuf;
	return 0;
}
