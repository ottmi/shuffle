#include "helper.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstdio>
#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;

int getMyId()
{
#if defined (_MPI)
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
#elif defined (_OPENMP)
    return omp_get_thread_num();
#else
    return 0;
#endif
}

int getNumOfCpus()
{
#if defined (_MPI)
    int num;
    MPI_Comm_size(MPI_COMM_WORLD, &num);
    return num;
#elif defined (_OPENMP)
    return omp_get_max_threads();
#else
    return 1;
#endif
}

double factorial(int n)
{
	double m = n;
	for (int i = n - 1; i > 1; i--)
		m *= i;

	return m;
}


string printTime(long t)
{
	stringstream s;
	if (t > 3600)
	{
		s << t / 3600 << ":" << setfill('0') << setw(2);
		t = t % 3600;
	}
	s << t / 60 << ":" << setfill('0') << setw(2) << t % 60;

	return s.str();
}

istream& safeGetline(istream& is, string& t)
{
	/* Courtesy of http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf */
	t.clear();
	istream::sentry se(is);
	streambuf* sb = is.rdbuf();

	for (;;)
	{
		int c = sb->sbumpc();
		switch (c)
		{
			case '\r':
				c = sb->sgetc();
				if (c == '\n')
					sb->sbumpc();
				return is;
			case '\n':
			case EOF:
				return is;
			default:
				t += (char) c;
		}
	}
}

string adjustString(string s, bool upercase)
{
	string r = "";

	for (unsigned int i = 0; i < s.length(); i++)
	{
		char c = s[i];
		if (c != '\t' && c != '\n' && c != '\r' && c != ' ')
		{
			if (upercase)
				r += toupper(c);
			else
				r += c;
		}
	}

	return (r);
}
