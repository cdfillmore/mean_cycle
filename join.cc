
//////////////////////////////////////////////////////
//                                                  //
//  Minimum Mean Weight Cycle Algorithm             //
//                                                  //
//  written by Christopher Fillmore                 //
//             Giridhar Shenoy                      //
//             Nadezhda Vassilyeva                  //
//                                                  //
//  to compile: g++ join.cc -o join -std=c++11      //
//  to execute: join w_test1.dmx                    //
//                                                  //
//////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <map>
#include <algorithm>

using namespace std;

//intersection of vectors
std::vector<int> intersection(std::vector<int> v1, std::vector<int> v2)
{
	std::vector<int> v3;
	
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	
	std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(v3));
	
	return v3;
}

//union of vectors
std::vector<int> union_1(std::vector<int> v1, std::vector<int> v2)
{
	std::vector<int> v3;
	
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	
	std::set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(v3));
	
	return v3;
}

//symmetric difference of vectors
std::vector<int> symmetric(std::vector<int> v1, std::vector<int> v2)
{
	std::vector<int> v3;
	
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	
	std::set_symmetric_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(v3));
	
	return v3;
}

//print a vector
void print_vector(std::vector<int> a)
{
	std::string out1="";
	for (int i=0; i<a.size(); i++)
	{
		out1.append(" ");
		out1.append(std::to_string(a.at(i)));
	}
	std::cout<<out1<<std::endl;
}

//read in the Graph (set the base index to 0 instead of 1)
std::map<int,std::map<int,float> > read_graph(char* file)
{
	FILE* fp;
	char line[4092];
	std::string summary;
	std::string tmp;
	int v_max;
	int e_max;
	int a;
	int b;
	float c;

	fp = fopen (file , "r");
	summary = fgets(line, sizeof(line), fp);
	std::stringstream ss;

	ss << summary;
	ss >> tmp >> tmp >> v_max >> e_max;
	
	std::map<int,std::map<int,float> > G;
	std::stringstream entry;

	for (int i=0; i<v_max; i++)
	{
		G[i] = std::map<int,float>({});
	}

	for (int i=0; i<e_max; i++)
	{
		fgets(line, sizeof(line), fp);
		ss << line;
		ss >> tmp >> a >> b >> c;
		G[a-1].insert(std::pair<int,float>(b-1,c));
		G[b-1].insert(std::pair<int,float>(a-1,c));
	}
	fclose(fp);
	return G;
}

//print the matching in dimacs format (fix base index to 1 again)
void output_dmx(std::vector<int> V)
{
	std::string tmp;
	int v_max = V.size();
	int e_max = 0;
	std::vector<int> visited = {};
	for (int i=0; i<v_max; i++)
	{
		if (std::find(visited.begin(),visited.end(),i) == visited.end())
		{
			e_max++;
			visited.push_back(V.at(i));
		}
	}
	std::cout<<"p edge "<<v_max<<e_max<<std::endl;
	for (int i=0; i<e_max; i++)
	{
		if (visited.at(i) > V.at(visited.at(i)))
		{
			std::cout<<"e "<<V.at(visited.at(i))+1<<" "<<visited.at(i)+1<<std::endl;
		}
		else
		{
			std::cout<<"e "<<visited.at(i)+1<<" "<<V.at(visited.at(i))+1<<std::endl;
		}
	}
	std::cout<<tmp<<std::endl;
}

void reduce(std::map<int,std::map<int,float>>& G)
{
	float m=G.begin()->second.begin()->second;
	for (auto const v : G)
	{
		for (auto const w : v.second)
		{
			if (m<w.second)
			{
				m = w.second;
			}
		}
	}
	for (auto v : G)
	{
		for (auto w : v.second)
		{
			G.at(v.first).at(w.first) = w.second - m;
		}
	}
}

//dijkstra's algorithm
std::pair<int,std::map<int,float>> dijk(std::map<int,std::map<int,float> > G, int s)
{
	std::map<int,float> l;
	std::vector<int> R = {};
	std::vector<int> V = {};
	int min;

	for (int i=0;i<G.size();i++)
	{
		V.push_back(i);
		l.insert(std::pair<int,float>(i,std::numeric_limits<float>::infinity()));
	}
	l[s] = 0;
	std::vector<int> iter = intersection(V,symmetric(R,V));
	while (iter.size() > 0 )
	{
		min = iter.at(0);
		for (int i=1;i<iter.size();i++)
		{
			if (l.at(iter.at(i)) < l.at(min))
			{
				min=iter.at(i);
			}
		}
		R = union_1(R,std::vector<int> {min});
		iter = intersection(V,symmetric(R,V));
		for (int i=0;i<iter.size();i++)
		{
			if (G.at(min).find(iter.at(i))!=G.at(min).end())
			{
				if (l.at(iter.at(i)) > (l.at(min) + G.at(min).at(iter.at(i))))
				{
					l.at(iter.at(i)) = l.at(min) + G.at(min).at(iter.at(i));
				}
			}
		}
	}
	return std::pair<int,std::map<int,float>> {s,l};
}

//metric closure
void metric_closure(std::map<int,std::map<int,float> > G)
{
	
}

//find minimum weight empty set join
void empty_tjoin(std::map<int,std::map<int,float> > G)
{
	
}

//minimum weight perfect matching
void MWPMP()
{

}

//main algorithm
void main_algorithm(std::map<int,std::map<int,float> > G, int mean)
{
	//initialize parameters
	int h=0;
	//reduce(G);
	std::cout<<dijk(G,0).second.at(1)<<std::endl;
}

//main function
int main(int argc, char* argv[])
{
	//check args for dmx file
	if (argc < 1)
	{
		cout << "no graph specified" << endl;
		return -1;
	}
	else
	{
		system("echo \"hi\"");
		//read Graph
		std::map<int,std::map<int,float> > G = read_graph(argv[1]);
		//set mean to zero
		int mean = 0;
		main_algorithm(G,mean);
		//run algorithm
		//std::vector<int> mu = find_max_mu(G,mean);
		//print the mu vector
		//std::cout<<mean<<std::endl;
		//output_dmx(mu);
		return 0;
	}
}

// The End //

