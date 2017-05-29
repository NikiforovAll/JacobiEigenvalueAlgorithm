#include "timer.h"
#include <iostream>

std::istream& operator>>(std::istream& is, timer** t)
{
	// Read timer tree from input file
	//
	// Format: N followed by N (name, elapsedtime) pairs then a list of N parents
	//

	int N;
	is >> N;

	// Read N lines of (name, elapsedtime) pairs
	timer* timers = new timer[N];
	for (int i = 0; i<N && !is.eof(); ++i)
	{
		std::string name;
		double value;
		is >> name >> value;

		timers[i].name = name;
		timers[i].elapsed = value;
	}

	// Build tree
	for (int i = 0; i<N; ++i)
	{
		int parent;
		is >> parent;

		if (parent != -1)
			timers[parent].add(&timers[i]);
	}

	*t = &timers[0];

	return is;
}

std::ostream& operator<<(std::ostream& os, timer& root)
{
	// Write timer tree to output file
	//
	// Format: N followed by N (name, elapsedtime) pairs then a list of N parents
	//
	typedef std::vector<std::pair<int, std::pair<std::string, double> > >
		vec_p_name_time;
	vec_p_name_time times;

	std::stack<std::pair<int, timer*> > Q;
	int lastparent = 0;

	Q.push(std::make_pair(-1, &root));

	while (Q.size())
	{
		timer* u = Q.top().second;
		times.push_back(std::make_pair(Q.top().first,
			std::make_pair(u->get_name(), u->get_elapsed())));
		Q.pop();

		// Add children to stack
		for (std::vector<timer*>::const_iterator it = u->get_children().begin();
			it != u->get_children().end();
			++it)
		{
			Q.push(std::make_pair(lastparent, *it));
		}
		lastparent++;
	}

	os << times.size() << std::endl;

	// Write timer data
	for (vec_p_name_time::iterator it = times.begin(); it != times.end(); ++it)
		os << it->second.first << " " << it->second.second << std::endl;

	// Write edges
	for (vec_p_name_time::iterator it = times.begin(); it != times.end(); ++it)
		os << it->first << std::endl;

	return os;
}


