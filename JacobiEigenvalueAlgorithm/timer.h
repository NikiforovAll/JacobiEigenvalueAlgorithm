#pragma once

#include <sys/time.h>
#include <sys/times.h>
#include <vector>
#include <string>
#include <stack>
#include <algorithm>
#include <fstream>
#include <iomanip>
//#ifdef enable_openmp
#include <omp.h>
//#endif

class timer {
public:
	timer()
		: calls(0), elapsed(0.0)
	{
	}
	timer(const std::string& name)
		: calls(0), elapsed(0.0), name(name)
	{
	}
	~timer()
	{
		// Each timer has sole ownership of its children so delete them.
		for (std::vector<timer*>::iterator it = children.begin();
			it != children.end();
			++it)
		{
			delete *it;
		}
	}

	void start()
	{
		++calls; begin = get_time();
	}

	void end()
	{
		elapsed += get_time() - begin;
	}

	void set_name(const std::string& n)
	{
		name = n;
	}

	void add(timer* p)
	{
		children.push_back(p);
	}

	const double& get_elapsed() const
	{
		return elapsed;
	}

	inline double get_time()
	{
		struct timeval  time_val;
		struct timezone time_zone;

		gettimeofday(&time_val, &time_zone);

		return (double)(time_val.tv_sec) + (double)(time_val.tv_usec) / 1000000.0;
	}

	double get(const std::string& n) const
	{
		if (name == n) return elapsed;
		for (std::vector<timer*>::const_iterator it = children.begin();
			it != children.end();
			++it)
		{
			double g = (*it)->get(n);
			if (g != -1.0)
				return g;
		}
		return -1.0;
	}

	const std::string& get_name() const { return name; }

	void print(std::ostream& os, int depth, double parentElapsed) const
	{
		// Do a depth first search printing child timers.

		// Indentation of name field and alignment of second column
		std::string indentation(depth * 2, ' ');
		int alignment = 28 - depth * 2 - name.size();
		std::string padding;
		for (int i = 0; i<alignment; ++i) padding += ' ';

		if (depth == 0)
			os << "name                        calls       elapsed     avg     % of parent\n"
			<< "-----------------------------------------------------------------------\n";

		os << indentation << name << padding
			<< std::fixed << std::setprecision(3)
			<< calls
			<< ((calls<1000) ? "\t" : "")	// padding
			<< "        " << elapsed
			<< "       " << (elapsed / static_cast<double>(calls))
			<< "    " << (elapsed / parentElapsed)*100.0
			<< '\n';

		for (std::vector<timer*>::const_iterator it = children.begin();
			it != children.end();
			++it)
		{
			(*it)->print(os, depth + 1, elapsed);
		}

	}

	//==========================================================================
	// Functions to write timer data to files 
	//==========================================================================
	void print_threads_v_elapsed() const
	{
		const std::string dir = "./timings/";
		const std::string postfix = "_elapsed.txt";
		std::string fname = dir + name + postfix;

		std::ofstream log(fname.c_str(), std::ios::app);
		#ifdef enable_openmp
		log << omp_get_max_threads() << '\t' << elapsed << '\n';
		#else
		log << 1 << '\t' << elapsed << '\n';
		#endif
		log.close();

		for (std::vector<timer*>::const_iterator it = children.begin();
			it != children.end();
			++it)
		{
			(*it)->print_threads_v_elapsed();
		}
	}

	#ifdef enable_openmp
	void print_threads_v_speedup(timer* serial) const
	{
		const std::string dir = "./timings/";
		const std::string postfix = "_speedup.txt";
		std::string fname = dir + name + postfix;

		std::ofstream log(fname.c_str(), std::ios::app);
		log << omp_get_max_threads() << '\t' << serial->get(name) / elapsed
			<< '\n';
		log.close();

		for (std::vector<timer*>::const_iterator it = children.begin();
			it != children.end();
			++it)
		{
			(*it)->print_threads_v_speedup(serial);
		}
	}

	/*	void print_threads_v_efficiency(timer* serial) const
	{
	const std::string dir = "./timings/";
	const std::string postfix = "_efficiency.txt";
	std::string fname = dir+name+postfix;

	std::ofstream log(fname.c_str(), std::ios::app);
	log << omp_get_max_threads() << '\t' << (serial->get(name)/elapsed)*(100.0/omp_get_max_threads())
	<< '\n';
	log.close();

	for(std::vector<timer*>::const_iterator it=children.begin();
	it != children.end();
	++it)
	{
	(*it)->print_threads_v_efficiency(serial);
	}
	}*/
#endif

	friend std::istream& operator>>(std::istream& is, timer** t);
	friend std::ostream& operator<<(std::ostream& os, timer& root);

private:
	const std::vector<timer*>& get_children() const { return children; }

	double begin;
	double elapsed;
	int calls;
	std::string name;
	std::vector<timer*> children;
};

std::istream& operator>>(std::istream& is, timer** t);
std::ostream& operator<<(std::ostream& os, timer& root);

