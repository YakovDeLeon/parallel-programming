#include "mpi.h"
#include <iostream>
#include <time.h>
#include <chrono>
#include <limits>
#include <string>

template<typename T>
void ShellSort(T* array, int n)
{
	T tmp;
	for (int step = n / 2; step > 0; step /= 2)
		for (int i = step; i < n; i++)
			for (int j = i - step; (j >= 0) && (array[j] > array[j + step]); j -= step)
			{
				tmp = array[j];
				array[j] = array[j + step];
				array[j + step] = tmp;
			}
}

template<typename T>
void FillingArray(T array[], int n, T min, T max) 
{
	srand(time(nullptr));
	for (int i = 0; i < n; i++)
	{
		T ri = (double)rand() / RAND_MAX;
		array[i] = min + (max - min)*ri;
	}
	std::cout << std::endl;
}


template <typename T>
void Comparator(T* arr, size_t size)
{
	for (size_t i = 1; i < size; i++)
	{
		if (arr[i] < arr[i - 1])
			std::swap(arr[i], arr[i - 1]);
	}
}

template <typename T>
void EvenSplitter(T* array, T* tmp1, T* tmp2, size_t size1, size_t size2)
{
	int a = 0;
	int b = 0;
	int i = 0;
	bool ChangeI = size1 % 2;

	while ((a < size1) && (b < size2))
	{

		if (tmp1[a] <= tmp2[b])
		{
			array[i] = tmp1[a];
			a += 2;
		}
		else
		{
			array[i] = tmp2[b];
			b += 2;
		}
		i += 2;
		if (i > size1 && ChangeI)
		{
			i -= 1;
			ChangeI = false;
		}
	}

	if (a == size1)
		for (int j = b; j < size2; j += 2, i += 2)
		{
			if (i > size1 && ChangeI)
			{
				i -= 1;
				ChangeI = false;
			}
			array[i] = tmp2[j];
		}
	else
		for (int j = a; j < size1; j += 2, i += 2)
		{
			if (i > size1 && ChangeI)
			{
				i -= 1;
				ChangeI = false;
			}
			array[i] = tmp1[j];
		}
}

template <typename T>
void OddSplitter(T* array, T* tmp1, T* tmp2, int size1, int size2)
{
	int a = 1;
	int b = 1;
	int i = 1;
	bool ChangeI = size1 % 2;

	while ((a < size1) && (b < size2))
	{

		if (tmp1[a] <= tmp2[b])
		{
			array[i] = tmp1[a];
			a += 2;
		}
		else
		{
			array[i] = tmp2[b];
			b += 2;
		}
		i += 2;
		if (i >= size1 && ChangeI)
		{
			i += 1;
			ChangeI = false;
		}
	}

	if (a == size1)
		for (int j = b; j < size2; j += 2, i += 2)
		{
			if (i >= size1 && ChangeI)
			{
				i += 1;
				ChangeI = false;
			}
			array[i] = tmp2[j];
		}
	else
		for (int j = a; j < size1; j += 2, i += 2)
		{
			if (i >= size1 && ChangeI)
			{
				i += 1;
				ChangeI = false;
			}
			array[i] = tmp1[j];
		}
}

template<typename T>
void ShowArray(T * pdArray, size_t unSizeArr, int ProcRank)
{
	for (size_t i = 0; i < unSizeArr; ++i)
	{
		std::cout << "ProcRank " << ProcRank << " has a folowing arr[" << i << "] = " << pdArray[i] << std::endl;
	}
}


template <typename T>
void BatcherMerge(T* arr1, T* arr2, T* res, size_t size1, size_t size2)
{
	memcpy(res, arr1, size1 * sizeof(T));
	memcpy(res + size1, arr2, size2 * sizeof(T));
	EvenSplitter(res, arr1, arr2, size1, size2);
	OddSplitter(res, arr1, arr2, size1, size2);
	Comparator(res, size1 + size2);
}

template <typename T>
void reallocate(T* &arr, size_t oldsize, size_t newsize)
{
	T* newarr = new T[newsize]{ 0 };
	memcpy(newarr, arr, oldsize * sizeof(T));
	delete[] arr;
	arr = newarr;
}

template <typename T>
int TREE_Betcher(T *sendbuf, T *recvbuf, int* counts, int sizearr, MPI_Datatype type, int root, MPI_Comm comm)
{
	MPI_Status status;
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int RankRecv, RankSend;
	T * ResBuf = new T[sizearr];
	int RecvArrSize = 0;
	int SendArrSize = counts[ProcRank];
	int newProcRank = (ProcRank - root + ProcNum) % ProcNum;
	int mask = 1;
	while (mask < ProcNum)
	{
		if ((newProcRank & mask) == 0)
		{
			RankSend = newProcRank | mask;
			if (RankSend < ProcNum)
			{
				RankSend = (RankSend + root) % ProcNum;
				MPI_Recv(&RecvArrSize, 1, MPI_INT, RankSend, 0, comm, &status);
				MPI_Recv(recvbuf, RecvArrSize, type, RankSend, 0, comm, &status);
				
				BatcherMerge(sendbuf, recvbuf, ResBuf, SendArrSize, RecvArrSize);
				reallocate(sendbuf, SendArrSize, SendArrSize + RecvArrSize);

				SendArrSize += RecvArrSize;
				memcpy(sendbuf, ResBuf, SendArrSize * sizeof(T));
				
			}
		}
		else
		{
			RankRecv = newProcRank&(~mask);
			RankRecv = (RankRecv + root) % ProcNum;
			MPI_Send(&SendArrSize, 1, MPI_INT, RankRecv, 0, comm);
			MPI_Send(sendbuf, SendArrSize, type, RankRecv, 0, comm);
			break;
		}
		mask = mask << 1;
	}
	if (ProcRank != root)
	{
		delete[] sendbuf;
		delete[] ResBuf;
	}

	if (ProcRank == root)
	{
		memcpy(recvbuf, ResBuf, sizearr * sizeof(T));
		delete[] sendbuf, ResBuf;
	}

	return 0;
}


int main(int argc, char *argv[])
{
	using namespace std;
	using namespace std::chrono;
	using time = chrono::steady_clock::time_point;
	int ProcNum, ProcRank;
	double *array = nullptr; 
	int size, root;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (argc == 5)
		root = atoi(argv[4]);
	else
		root = 0;
	if (ProcRank == root)
	{
		if (argc == 5)
			size = atoi(argv[1]);
		else size = 10;
	}
	MPI_Bcast(&size, 1, MPI_INT, root, MPI_COMM_WORLD);

	high_resolution_clock::time_point startTime;
	high_resolution_clock::time_point endTime;
	
	if (ProcRank == root)
	{
		array = new double[size];
		FillingArray(array, size, atoi(argv[2]), atoi(argv[3]));
		//ShowArray(array, size);
		
	}
	
	int * sendcounts = new int[ProcNum];
	int * displs = new int[ProcNum];
	if (ProcRank == root)
	{
		int nSendSum = 0;
		int nRecvSum = 0;
		int nRemSizePerProc = size % ProcNum;
		for (size_t i = 0; i < ProcNum; i++)
		{
			sendcounts[i] = size / ProcNum;
			if (nRemSizePerProc > 0)
			{
				sendcounts[i]++;
				nRemSizePerProc--;
			}
			displs[i] = nRecvSum;
			nRecvSum += sendcounts[i];
		}
		startTime = chrono::high_resolution_clock::now();
	}

	MPI_Bcast(sendcounts, ProcNum, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(displs, ProcNum, MPI_INT, root, MPI_COMM_WORLD);

	double *recvbuf = new double[sendcounts[ProcRank]];

	MPI_Scatterv(array, sendcounts, displs, MPI_DOUBLE, recvbuf, sendcounts[ProcRank], MPI_DOUBLE, root, MPI_COMM_WORLD);

	ShellSort(recvbuf, sendcounts[ProcRank]);

	double * res = new double[size];
	TREE_Betcher(recvbuf, res, sendcounts, size, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
	if (ProcRank == root)
	{
		//std::cout << "Result for parallel version: " << std::endl;
		//ShowArray(res, size);
		endTime = chrono::high_resolution_clock::now();
		auto WorkTimeParallel = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
		cout << "Time for parallel version: " << WorkTimeParallel << "ns" << endl;
		startTime = chrono::high_resolution_clock::now();
		ShellSort(array, size);
		//std::cout << "Result for linear version: " << std::endl;
		//ShowArray(array, size);
		endTime = chrono::high_resolution_clock::now();
		auto WorkTimeSerial = chrono::duration_cast<chrono::nanoseconds>(endTime - startTime).count();
		cout << "Time for linear version: " << WorkTimeSerial << "ns" << endl;
		string compare = (WorkTimeParallel < WorkTimeSerial ? "faster" : "slower");
		cout << "Parallel version " << compare << " linear version!" << std::endl;
		bool flag;
		for (int i = 0; i < size; i++)
		{
			flag = true;
			if (res[i] != array[i])
			{
				std::cout << "The parallel and linear version not matching" << std::endl;
				flag = false;
				break;
			}
		}
		
		if (flag)
		{
			std::cout << "The parallel and linear version matching" << std::endl;
		}
	}
	MPI_Finalize();
	delete[] array, sendcounts, displs, recvbuf, res;
	return 0;
}
